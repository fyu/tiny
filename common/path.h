#ifndef FURRY_COMMON_PATH
#define FURRY_COMMON_PATH

#include <string>
#include <vector>

namespace furry {

bool MakeDir(const std::string &dir);
bool mkdir(const std::string &dir);

bool remove_file(const std::string &filename);

// std::string SplitFilename(const std::string &path);

bool Exists(const std::string &path);
bool exists(const std::string &path);

std::string GetExt(const std::string &filename);

bool ListDir(const std::string &dirname, const std::string &ext,
             std::vector<std::string> *filenames);

std::string GetExt(const std::string &filename);

std::string get_ext(const std::string &filename);
std::string get_stem(const std::string &filename);
std::string get_dir(const std::string &filename);
std::string add_suffix(const std::string &filename, const std::string &suffix);
std::string add_index(const std::string &filename, int index,
                      const std::string &pattern = "%02d");

} // furry

#endif // FURRY_COMMON_PATH
