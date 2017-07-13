#include "furry/common/protobuf.h"

#include <fcntl.h>

#include <glog/logging.h>

namespace furry {

bool ReadTxtFileToProto(const std::string& filename,
                     google::protobuf::Message* proto) {
  int file = open(filename.c_str(), O_RDONLY);
  if (file == -1) {
    LOG(ERROR) << "Failed to open " << filename << " to read";
    return false;
  }
  google::protobuf::io::FileInputStream file_stream(file);
  return google::protobuf::TextFormat::Parse(&file_stream, proto);
}

bool WriteProtoToTxtFile(const std::string &filename,
                         const google::protobuf::Message &proto) {
  int file = open(filename.c_str(), O_WRONLY | O_CREAT,
                  S_IRWXU | S_IRWXG | S_IRWXO);
  if (file == -1) {
    LOG(ERROR) << "Failed to open " << filename << " to write";
    return false;
  }
  google::protobuf::io::FileOutputStream file_stream(file);
  return google::protobuf::TextFormat::Print(proto, &file_stream);
}

bool read_txt_proto(const std::string &filename,
                    google::protobuf::Message* proto) {
  return ReadTxtFileToProto(filename, proto);
}
bool write_txt_proto(const std::string &filename,
                     const google::protobuf::Message &proto) {
  return WriteProtoToTxtFile(filename, proto);
}

} // furry
