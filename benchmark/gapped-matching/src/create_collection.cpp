#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_mapper.hpp>
#include <iostream>

#include "utils.hpp"
#include "constants.hpp"
#include "collection.hpp"
#include "logging.hpp"

typedef struct cmdargs {
    std::string input_file;
    std::string collection_dir;
} cmdargs_t;

void print_usage(const char* program)
{
    fprintf(stdout, "%s -i <input file> -c <col dir>\n", program);
    fprintf(stdout, "where\n");
    fprintf(stdout, "  -i <input file>      : the input file.\n");
    fprintf(stdout, "  -c <collection dir>  : the collection dir.\n");
};

cmdargs_t parse_args(int argc, const char* argv[])
{
    cmdargs_t args;
    int op;
    args.input_file = "";
    args.collection_dir = "";
    while ((op = getopt(argc, (char* const*)argv, "i:c:")) != -1) {
        switch (op) {
            case 'i':
                args.input_file = optarg;
                break;
            case 'c':
                args.collection_dir = optarg;
                break;
        }
    }
    if (args.collection_dir == "" || args.input_file == "") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

int main(int argc, const char* argv[])
{
    log::start_log(argc, argv);
    cmdargs_t args = parse_args(argc, argv);

    /* (0) check if input file exists */
    std::ifstream ifs(args.input_file);
    if (!ifs) {
        LOG(FATAL) << "Error opening input file '" << args.input_file << "'";
    }

    /* (1) create collection dir */
    LOG(INFO) << "Creating collection directory structure in '" << args.collection_dir << "'";
    utils::create_directory(args.collection_dir);

    /* (2) create sub dirs */
    utils::create_directory(args.collection_dir+"/tmp");
    utils::create_directory(args.collection_dir+"/index");
    utils::create_directory(args.collection_dir+"/patterns");
    utils::create_directory(args.collection_dir+"/results");

    /* (3) copy file to sdsl format */
    auto buf = sdsl::write_out_buffer<0>::create(args.collection_dir+"/"+ consts::KEY_PREFIX + consts::KEY_TEXT);
    {
        gm_timer tm("COPY TO SDSL FORMAT",true);
        std::copy(std::istream_iterator<uint8_t>(ifs),
                  std::istream_iterator<uint8_t>(),
                  std::back_inserter(buf));
    }
    {
        gm_timer tm("BIT COMPRESS",true);
        sdsl::util::bit_compress(buf);
    }

    LOG(INFO) << "num_syms = " << buf.size();
    LOG(INFO) << "log2(sigma) = " << (int) buf.width();
    LOG(INFO) << "text size = " << (buf.width()*buf.size())/(8*1024*1024) << " MiB";

    return EXIT_SUCCESS;
}