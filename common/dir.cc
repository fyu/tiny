#include "dir.h"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <pwd.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>

namespace ba = boost::algorithm;

namespace furry
{
std::string
get_home_dir()
{
  struct passwd *pws;
  pws = getpwuid(geteuid());
  return pws->pw_dir;
}

std::string HomeDir() {
  struct passwd *pws;
  pws = getpwuid(geteuid());
  return pws->pw_dir;
}

#ifdef _WIN32
static const char *separator = "\\";
#else
static const char *separator = "/";
#endif

std::vector<std::string>
split_ext(const std::string & file_name)
{
  std::vector<std::string> result(2);
  auto pos = file_name.rfind('.');
  result[0] = file_name.substr(0, pos);
  if (pos != std::string::npos)
    result[1] = file_name.substr(pos+1);

  return result;
}

std::vector<std::string>
split_base(const std::string & file_name)
{
  std::vector<std::string> result(2);
  auto pos = file_name.find_last_of("\\/");
  if (pos == std::string::npos)
  {
    result[1] = file_name;
  }
  else
  {
    result[0] = file_name.substr(0, pos+1);
    result[1] = file_name.substr(pos+1);
  }
  return result;
}

std::string
base_name(const std::string & file_name)
{
  auto pos = file_name.find_last_of("\\/");
  if (pos != std::string::npos)
    return file_name.substr(pos+1);
  else
    return file_name;
}

std::vector<std::string>
ListDir(const std::string & dir_name,
          const std::string & ext)
{
  std::vector<std::string> names;
  DIR *pDIR;
	struct dirent *entry;
  pDIR = opendir(dir_name.c_str());
  if(pDIR){
		while((entry = readdir(pDIR))){
			if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 &&
          (ext == "" || split_ext(entry->d_name)[1] == ext))
      {
        names.push_back(entry->d_name);
      }
		}
    closedir(pDIR);
	}

  return names;
}

std::vector<std::string>
SplitPath(const std::string &path) {
  std::vector<std::string> path_parts;
  ba::split(path_parts, path,
            ba::is_any_of(separator),
            ba::token_compress_on);
  return path_parts;
}

std::string
PathToId(const std::string &path) {
  std::list<std::string> path_parts;
  ba::split(path_parts, path, ba::is_any_of("./"), ba::token_compress_on);
  if (path_parts.size() != 0 && path_parts.back().size() == 0)
    path_parts.pop_back();
  if (path_parts.size() != 0 && path_parts.front().size() == 0)
    path_parts.pop_front();
  if (path_parts.size() == 0)
    return "";
  else {
    std::stringstream ss;
    std::copy(path_parts.begin(),
              path_parts.end(),
              std::ostream_iterator<std::string>(ss, "_"));
    auto id = ss.str();
    return id.substr(0, id.size() - 1);
  }
}

std::string GetExtension(const std::string &name) {
  size_t pos = name.rfind('.');
  if (pos == std::string::npos) return "";
  else return name.substr(pos);
}

std::string AddExtension(const std::string &name, const std::string &ext) {
  boost::filesystem::path p(name);
  if (p.extension() != ext)
    p += ext;
  return p.string();
}

// Boost version
// void MakeDir(const std::string &dir) {
//   std::vector<std::string> pieces;
//   boost::split(pieces, dir, boost::is_any_of("/"), boost::token_compress_on);
//   boost::filesystem::path p;
//   if (pieces.size() > 0 && dir[0] == '/')
//     p /= '/' + pieces[0];
//   else
//     p /= pieces[0];
//   boost::filesystem::create_directory(p);
//   for (size_t i = 1; i < pieces.size(); ++i) {
//     p /= pieces[i];
//     boost::filesystem::create_directory(p);
//   }
// }

} // furry
