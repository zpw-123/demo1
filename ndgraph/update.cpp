#define ENABLE_LOCK 0
//#define WEIGHTED 0
#define VERIFY 0
#define OPENMP 1
#include <assert.h>
#include <stdlib.h>


#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <vector>

#include "graph.hpp"
#include "io_util.h"
#include "parallel.h"
#include "rmat_util.h"
#include "sys/time.h"


#include "graph_converter.hpp"

std::vector<uint32_t> get_random_permutation_new(uint32_t num) {
  std::vector<uint32_t> perm(num);
  std::vector<uint32_t> vec(num);

  for (uint32_t i = 0; i < num; i++) vec[i] = i;

  uint32_t cnt{0};
  while (vec.size()) {
    uint32_t n = vec.size();
    srand(time(NULL));
    uint32_t idx = rand() % n;
    uint32_t val = vec[idx];
    std::swap(vec[idx], vec[n - 1]);
    vec.pop_back();
    perm[cnt++] = val;
  }
  return perm;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    printf("Invalid arguments.\n");
    return -1;
  }

  string file_name = argv[1];
  //string file_name = "/home/zpw/HME-Quartz-broadwell-master/dataset/soc-LiveJournal1.txt";
  std::string src, dest;
  // read edges as source and destination
  // int cnt = 0;
  struct timeval start, end;
  struct timezone tzp;

  // initialize graph
  //uint32_t num_nodes=4847571;
  //uint64_t num_edges=68993760;
  //uint32_t num_nodes=9;
  //uint64_t num_edges=30;
  bool wei1 = false;
  bool sorted   = true;
  bool skip     = false;
  size_t blocksize = 256 * 1024 * 1024;
  
  auto static_query_blocksize = [blocksize](vid_t nvertices) { return blocksize; };
  std::function<size_t(vid_t nvertices)> query_blocksize;
  query_blocksize = static_query_blocksize;
  graph_converter converter(remove_extension(argv[1]), wei1, sorted);
  convert(file_name, converter, query_blocksize, skip);
  
  std::string base_name = remove_extension(argv[1]);

    /* graph meta info */
  unsigned int num_nodes;
  unsigned long int num_edges;
  load_graph_meta(base_name, &num_nodes, &num_edges, wei1);
  std::cout<<num_nodes<<" "<<num_edges<<std::endl;

  std::string vert_block_name = get_vert_blocks_name(base_name, blocksize);
  std::string edge_block_name = get_edge_blocks_name(base_name, blocksize);
  std::string weight_name = get_weights_name(base_name);
  std::string beg_pos_name = get_beg_pos_name(base_name);
  std::string csr_name = get_csr_name(base_name);
  std::cout<<vert_block_name<<std::endl;


  std::vector<vid_t> vblocks = load_graph_blocks<vid_t>(vert_block_name);
  std::vector<eid_t> eblocks = load_graph_blocks<eid_t>(edge_block_name);

  int vertdesc = open(beg_pos_name.c_str(), O_RDONLY);
  int edgedesc = open(csr_name.c_str(), O_RDWR);

  pre_block_t block;
  bid_t nblocks = vblocks.size() - 1;
  std::cout<<nblocks<<std::endl;

  pair_uint *edges = (pair_uint *)malloc((num_edges) * sizeof(pair_uint));
  std::cout<<"edges init"<<std::endl;
  unsigned long int start1=0;
  unsigned int a;
  unsigned int b;
  for (bid_t blk = 0; blk < nblocks; blk++)
  {

    block.nverts = vblocks[blk + 1] - vblocks[blk];
    block.nedges = eblocks[blk + 1] - eblocks[blk];
    block.start_vert = vblocks[blk];
    block.start_edge = eblocks[blk];

    block.beg_pos = (eid_t *)realloc(block.beg_pos, (block.nverts + 1) * sizeof(eid_t));
    block.csr = (vid_t *)realloc(block.csr, block.nedges * sizeof(vid_t));
    std::cout<<"beg_pos init"<<std::endl;

    load_block_range(vertdesc, block.beg_pos, (block.nverts + 1), block.start_vert * sizeof(eid_t));
    load_block_range(edgedesc, block.csr, block.nedges, block.start_edge * sizeof(vid_t));
    std::cout<<"load_block_range"<<std::endl;

    for(uint32_t ii=0;ii<block.nverts;ii++)
    {
      //由于是累加的
      //std::cout<<ii<<std::endl;
      for(uint64_t jj=block.beg_pos[ii];jj<block.beg_pos[ii+1];jj++)
      {
        if(blk>0)
        {
          a=blk*(vblocks[1] - vblocks[0])+ii;
        }
        else{
          a=ii;
        }
        
        b=block.csr[jj%block.nedges];
        edges[jj] = {a,b};
        if(jj>=num_edges)break;
        
      }
    }

  }



  std::cout<<"edges okk"<<std::endl;

  //pair_uint *edges =
  //    get_edges_from_file(file_name.c_str(), &num_edges, &num_nodes);

  

  Graph graph(num_nodes);
  std::vector<uint32_t> new_srcs(num_edges);
  std::vector<uint32_t> new_dests(num_edges);
  for (uint32_t i = 0; i < num_edges; i++) {
    new_srcs[i] = edges[i].x;
    new_dests[i] = edges[i].y;
    //cout<<new_srcs[i]<<" "<<new_dests[i]<<endl;
  }
  auto perm = get_random_permutation_new(num_edges);  //
  
  std::cout << "Insert edges" << std::endl;
  gettimeofday(&start, &tzp);
  ofstream outfile;
  outfile.open("data.txt", ios::binary | ios::app | ios::in | ios::out);
  outfile<<tzp.tz_dsttime<<"  minute:   "<<tzp.tz_minuteswest<<"\n";
  outfile<<"miao: "<<start.tv_sec<<"  umiao:   "<<start.tv_usec<<"\n";
  std::cout << start.tv_sec << std::endl;
  std::cout << start.tv_usec << std::endl;
  std::cout<<tzp.tz_minuteswest<<std::endl;
  std::cout<<tzp.tz_dsttime<<std::endl;
  graph.add_edge_batch(new_srcs.data(), new_dests.data(), num_edges, perm);
  gettimeofday(&end, &tzp);
  outfile<<tzp.tz_dsttime<<"  minute:   "<<tzp.tz_minuteswest<<"\n";
  outfile<<"miao: "<<end.tv_sec<<"  umiao:   "<<end.tv_usec<<"\n";
  std::cout << end.tv_sec << std::endl;
  std::cout << end.tv_usec << std::endl;
  std::cout<<tzp.tz_minuteswest<<std::endl;
  std::cout<<tzp.tz_dsttime<<std::endl;
  
}