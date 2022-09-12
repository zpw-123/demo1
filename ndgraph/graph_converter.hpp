#ifndef _GRAPH_CONVERTER_H_
#define _GRAPH_CONVERTER_H_

#include <string>
#include <vector>
#include <fstream>
#include <functional>
#include "graph_buffer.hpp"

#include "constants.hpp"
#include "types.hpp"
#include "logger.hpp"
#include "util.hpp"
#include "io.hpp"
#include "precompute.hpp"
#include "split.hpp"


/** This file defines the data structure that contribute to convert the text format graph to some specific format */

class base_converter {
public:
    base_converter() = default;
    ~base_converter() { }
    virtual void initialize() { };
    virtual void convert(vid_t from, vid_t to) { };
    virtual void finalize() { };
};


/**
 * Convert the text format graph dataset into the csr format. the dataset may be very large, so we may split the csr data
 * into multiple files.
 * `fnum` records the number of files used to store the csr data.
 * `beg_pos` : records the vertex points the position of csr array
 * `csr` : records the edge destination
 * `deg` : records the degree of each vertex
 *
 * `adj` : records the current vertex neighbors
 * `curr_vert` : the current vertex id
 * `max_vert` : records the max vertex id
 * `rd_edges` : the number of edges that has been read, use to split into multiply files.
 */
class graph_converter : public base_converter {
private:
    int fnum;
    bool _weighted;
    bool _sorted;
    graph_buffer<eid_t> beg_pos;
    graph_buffer<vid_t> csr;
    graph_buffer<vid_t> deg;
    graph_buffer<real_t> weights;

    std::vector<vid_t> adj;
    std::vector<real_t> adj_weights;
    vid_t curr_vert;
    vid_t max_vert;

    /* for buffer */
    vid_t buf_vstart; /* current buffer start vertex */
    eid_t buf_estart; /* current buffer start edge */
    eid_t rd_edges;
    eid_t csr_pos;

    std::string base_name;

    void setup_output(const std::string& input) {
        base_name = input;
    }

    void setup_output(const std::string& path, const std::string& dataset) {
        base_name = path + dataset;
    }

    void flush_beg_pos() {
        std::string name = get_beg_pos_name(base_name);
        appendfile(name, beg_pos.buffer_begin(), beg_pos.size());
        beg_pos.clear();
    }

    void flush_csr() {
        std::string name = get_csr_name(base_name);
        appendfile(name, csr.buffer_begin(), csr.size());
        rd_edges += csr.size();
        buf_estart += csr.size();
        csr.clear();
    }

    void flush_weights() {
        std::string name = get_weights_name(base_name);
        appendfile(name, weights.buffer_begin(), weights.size());
        weights.clear();
    }

    void sync_buffer() {
        for(auto & dst : adj) csr.push_back(dst);
        if(_weighted) {
            for(const auto & w : adj_weights) weights.push_back(w);
        }
        csr_pos += adj.size();
        beg_pos.push_back(csr_pos);
        deg.push_back(adj.size());
    }

    void sync_zeros(vid_t zeronodes) {
        while(zeronodes--){
            if(beg_pos.full()) flush_buffer();
            beg_pos.push_back(csr_pos);
            deg.push_back(0);
        }
    }
public:
    graph_converter() = delete;
    graph_converter(const std::string& path, bool weighted = false, bool sorted = false) {
        fnum = 0;
        beg_pos.alloc(VERT_SIZE);
        csr.alloc(EDGE_SIZE);
        deg.alloc(VERT_SIZE);
        curr_vert = max_vert = buf_vstart = buf_estart = rd_edges = csr_pos = 0;
        setup_output(path);
        _weighted = weighted;
        _sorted = sorted;
        if(_weighted) {
            weights.alloc(EDGE_SIZE);
        }
    }
    graph_converter(const std::string& folder, const std::string& dataset, bool weighted = false, bool sorted = false) {
        fnum = 0;
        beg_pos.alloc(VERT_SIZE);
        csr.alloc(EDGE_SIZE);
        deg.alloc(VERT_SIZE);
        curr_vert = max_vert = buf_vstart = buf_estart = rd_edges = csr_pos = 0;
        setup_output(folder, dataset);
        _weighted = weighted;
        _sorted = sorted;
        if(_weighted) {
            weights.alloc(EDGE_SIZE);
        }
    }
    graph_converter(const std::string& path, size_t vert_size, size_t edge_size, bool weighted = false, bool sorted = false) {
        fnum = 0;
        beg_pos.alloc(vert_size);
        csr.alloc(edge_size);
        deg.alloc(vert_size);
        curr_vert = max_vert = buf_vstart = buf_estart = rd_edges = csr_pos = 0;
        setup_output(path);
        _weighted = weighted;
        _sorted = sorted;
        if(_weighted) {
            weights.alloc(EDGE_SIZE);
        }
    }
    ~graph_converter() {
        beg_pos.destroy();
        csr.destroy();
        deg.destroy();
    }

    void initialize() {
        beg_pos.clear();
        csr.clear();
        deg.clear();
        curr_vert = max_vert = buf_vstart = buf_estart = rd_edges = csr_pos = 0;
        beg_pos.push_back(0);
    }

    void convert(vid_t from, vid_t to, real_t *weight) {
        max_vert = max_value(max_vert, from);
        max_vert = max_value(max_vert, to);

        if(from == curr_vert) {
            adj.push_back(to);
            if(_weighted) adj_weights.push_back(*weight);
        }
        else {
            if(csr.test_overflow(adj.size()) || beg_pos.full() ) {
                flush_buffer();
            }
            if(adj.size() > EDGE_SIZE) {
                logstream(LOG_ERROR) << "Too small memory capacity with EDGE_SIZE = " << EDGE_SIZE << " to support larger out degree = " << adj.size() << std::endl;
                assert(false);
            }
            sync_buffer();
            if(from - curr_vert > 1) sync_zeros(from - curr_vert - 1);
            curr_vert = from;
            adj.clear();
            adj.push_back(to);
            if(_weighted) {
                adj_weights.clear();
                adj_weights.push_back(*weight);
            }
        }
    }

    void flush_buffer() {
        logstream(LOG_INFO) << "Buffer : [ " << buf_vstart << ", " <<  buf_vstart + deg.size() << " ), csr position : [ " << buf_estart << ", " << csr_pos << " )" << std::endl;
        if(!csr.empty()) flush_csr();
        if(!beg_pos.empty()) flush_beg_pos();
        if(_weighted && !weights.empty()) flush_weights();
    }

    void finalize() {
        sync_buffer();
        if(max_vert > curr_vert) sync_zeros(max_vert - curr_vert);
        flush_buffer();

        /** write the graph meta data into meta file */
        std::string metafile = get_meta_name(base_name);
        auto metastream = std::fstream(metafile.c_str(), std::ios::out | std::ios::binary);
        vid_t nvertices = max_vert + 1;
        eid_t nedges = csr_pos;
        metastream.write((char *)&nvertices, sizeof(vid_t));
        metastream.write((char *)&nedges, sizeof(eid_t));
        metastream.close();

        logstream(LOG_INFO) << "nvertices = " << max_vert + 1 << ", nedges = " << csr_pos << ", files : " << fnum + 1 << std::endl;
    }

    int get_fnum() { return this->fnum + 1; }
    std::string get_output_filename() const { return base_name; }
    bool is_weighted() const { return _weighted; }

    vid_t get_nvertices() { return this->max_vert + 1; }

    bool need_sorted() const { return _sorted; }
};

void convert(std::string filename, graph_converter &converter, std::function<size_t(vid_t nvertices)> query_blocksize, bool skip = false) {

    std::string base_name = remove_extension(filename);
    bool reprocessed = true;
    if(test_dataset_processed_exists(base_name)) {
        if(skip) {
            reprocessed = false;
        } else {
            logstream(LOG_INFO) << "Find preprocessed data, do you want to reprocessed(y/n)?" << std::endl;
            std::string val;
            std::cin >> val;
            if (!val.empty() && (val[0] == 'n' || val[0] == 'N')) reprocessed = false;
        }
    }

    if (reprocessed) {
        delete_processed_dataset(base_name);
        FILE *fp = fopen(filename.c_str(), "r");
        assert(fp != NULL);

        logstream(LOG_INFO) << "start to convert the " << filename << std::endl;
        converter.initialize();
        char line[1024];
        size_t rdlines = 0;
        while (fgets(line, 1024, fp) != NULL)
        {
            if (line[0] == '#')
                continue;
            if (line[0] == '%')
                continue;

            char *t1, *t2, *t3;
            rdlines++;
            t1 = strtok(line, "\t, ");
            t2 = strtok(NULL, "\t, ");
            t3 = strtok(NULL, "\t, ");
            if (t1 == NULL || t2 == NULL)
            {
                logstream(LOG_ERROR) << "Input file is not the right format. Expected <from> <to>" << std::endl;
                assert(false);
            }
            vid_t from = atoi(t1);
            vid_t to = atoi(t2);
            if (from == to)
                continue;
            if (converter.is_weighted())
            {
                assert(t3 != NULL);
                real_t w = static_cast<real_t>(atof(t3));
                converter.convert(from, to, &w);
            }
            else
            {
                converter.convert(from, to, NULL);
            }
            if(rdlines % 1000000 == 0) {
                logstream(LOG_INFO) << "readlines : " << rdlines << std::endl;
            }
        }
        converter.finalize();
        logstream(LOG_INFO) << "finish to convert the " << filename << std::endl;
    }

    vid_t nvertices;
    eid_t nedges;
    load_graph_meta(base_name, &nvertices, &nedges, false);
    size_t blocksize = query_blocksize(nvertices);
    
    std::string folder = get_dataset_block_folder(base_name, blocksize);
    if(!test_folder_exists(folder)) {
        sowalker_mkdir(folder.c_str());
    }

    bool block_data_exist = test_dataset_block_data_exists(base_name,  blocksize);
    bool regenerate = true;
    if(!reprocessed && block_data_exist) {
        if(skip) {
            regenerate = false;
        } else {
            logstream(LOG_INFO) << "Find preprocessed block data (edge.block, vert.block, exp), do you want to regenerate(y/n)?" << std::endl;
            std::string val;
            std::cin >> val;
            if (!val.empty() && (val[0] == 'n' || val[0] == 'N')) regenerate = false;
        }
    }

    if(reprocessed || regenerate) {
        delete_processed_block_data(base_name, blocksize);
        /* split the data into multiple blocks */
        split_blocks(base_name, 0, blocksize);

        if(converter.need_sorted()) {
            sort_vertex_neighbors(base_name, blocksize, converter.is_weighted());
        }

        /* make the expected walk length */
        calc_expected_walk_length(base_name, blocksize, 10);
    }
}

#endif
