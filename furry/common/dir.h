#ifndef FURRY_COMMON_DIR
#define FURRY_COMMON_DIR

#include <string>
#include <vector>
#include <dirent.h>

#include <boost/filesystem.hpp>

#include "furry/common/cast.h"

namespace furry
{
std::string
get_home_dir();

std::string HomeDir();

std::vector<std::string>
split_ext(const std::string& file_name);

std::vector<std::string>
split_base(const std::string& file_name);

std::string
base_name(const std::string& file_name);

std::vector<std::string>
SplitPath(const std::string &path);

// Convert a path to an identifier to this path
std::string
PathToId(const std::string &path);

std::vector<std::string>
ListDir(const std::string& dir_name,
         const std::string& ext = "");

std::string AddExtension(const std::string &name, const std::string &ext);

template <typename Result, typename... Args> Result
SaveInPath(std::function<Result (const std::string&)> save_action,
           boost::filesystem::path current_path,
           Args... path);

template <typename Result, typename T> inline Result
SaveInPath(std::function<Result (const std::string&)> save_action,
           boost::filesystem::path current_path,
           T file_name) {
  current_path /= To<std::string>(file_name);
  return save_action(current_path.string());
}

template <typename Result, typename T, typename... Args> inline
Result
SaveInPath(std::function<Result (const std::string&)> save_action,
           boost::filesystem::path current_path,
           T first_level,
           Args... rest) {
  current_path /= To<std::string>(first_level);
  boost::filesystem::create_directory(current_path);
  return SaveInPath(save_action, current_path, rest...);
}

template <typename Result, typename... Args> Result
SaveInPath(std::function<Result (const std::string&)> save_action,
           Args... path) {
  return SaveInPath(save_action,
                    boost::filesystem::path(),
                    path...);
}

template <typename... Args> void
SaveInPath(std::function<void (const std::string&)> save_action,
           Args... path) {
  SaveInPath<void>(save_action, path...);
}

} // furry

// Non-boost functions are in this header
#include "furry/common/path.h"

#endif // FURRY_COMMON_PATH
