#ifndef FURRY_COMMON_PROTOBUF_H
#define FURRY_COMMON_PROTOBUF_H

#include <string>

#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

namespace furry {

bool ReadTxtFileToProto(const std::string &filename,
                        google::protobuf::Message* proto);

bool WriteProtoToTxtFile(const std::string &filename,
                         const google::protobuf::Message &proto);

bool read_txt_proto(const std::string &filename,
                    google::protobuf::Message* proto);
bool write_txt_proto(const std::string &filename,
                     const google::protobuf::Message &proto);


} // furry

#endif // FURRY_COMMON_PROTOBUF_H
