#include <iostream>
#include "utils.h"

void block_convert(std::ifstream& reader, std::ofstream& writer,
                   _s32* read_buf, _s32* write_buf, _u64 npts, _u64 ndims) {
  reader.read((char*) read_buf,
              npts * (ndims * sizeof(_s32)));
  for (_u64 i = 0; i < npts; i++) {
    memcpy(write_buf+i*(ndims+1), &ndims, sizeof(unsigned));
    memcpy(write_buf + i*(ndims+1)+1, (read_buf + i * (ndims)),
           ndims * sizeof(_s32));
  }
  writer.write((char*) write_buf, npts * (ndims *sizeof(_s32) + sizeof(unsigned)));
}


int main(int argc, char** argv){
    if (argc != 3) {
        std::cout << argv[0] << " input_ibin output_ivecs" << std::endl;
        exit(-1);
    }
    std::ifstream reader(argv[1], std::ios::binary);
    int           npts_s32;
    int           ndims_s32;
    reader.read((char*) &npts_s32, sizeof(_s32));
    reader.read((char*) &ndims_s32, sizeof(_s32));
    size_t npts = npts_s32;
    size_t ndims = ndims_s32;
    std::cout << "Dataset: #pts = " << npts << ", # dims = " << ndims
            << std::endl;

    _u64 blk_size = 131072;
    _u64 nblks = ROUND_UP(npts, blk_size) / blk_size;
    std::cout << "# blks: " << nblks << std::endl;
    std::ofstream writer(argv[2], std::ios::binary);
    _s32* read_buf = new _s32[npts * (ndims)];
    _s32* write_buf = new _s32[npts * (ndims+1)];

    for (_u64 i = 0; i < nblks; i++) {
        _u64 cblk_size = std::min(npts - i * blk_size, blk_size);
        block_convert(reader, writer, read_buf, write_buf, cblk_size, ndims);
        std::cout << "Block #" << i << " written" << std::endl;
    }

    delete[] read_buf;
    delete[] write_buf;

    reader.close();
    writer.close();
}