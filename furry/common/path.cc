#include "furry/common/path.h"

#include <cstdlib>
#include <regex>
#include <string.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <glog/logging.h>

#include "furry/common/str.h"

namespace {
// Not supported yet
// const std::regex unix_file_name_exp = std::regex("^[A-Za-z0-9\\-\\._]+$");
}

namespace furry {

bool MakeDir(const std::string &dir) {
  // LOG(INFO) << "mkdir " << dir;
  std::string cmd = "mkdir -m 774 -v -p " + dir;
  return system(cmd.c_str()) == 0;
}

bool mkdir(const std::string &dir) {
  return MakeDir(dir);
}

bool remove_file(const std::string &filename) {
  std::string cmd = "rm -f " + filename;
  return system(cmd.c_str()) == 0;
}

bool Exists(const std::string &path) {
  struct stat st;
  if (stat(path.c_str(), &st) == 0)
    return true;
  else
    return false;
}

bool exists(const std::string &path) {
  return Exists(path);
}

std::string GetExt(const std::string &filename) {
  auto p = filename.rfind('.');
  if (p == std::string::npos) return "";
  std::string ext = filename.substr(p);
  if (ext.find_first_of("\\/") == std::string::npos) return ext;
  else return "";
}

std::string get_ext(const std::string &filename) {
  return GetExt(filename);
}

std::string get_stem(const std::string &filename) {
  auto p = filename.rfind('.');
  if (p == std::string::npos) return filename;
  return filename.substr(0, p);
}

std::string get_dir(const std::string &filename) {
  auto p = filename.rfind('/');
  if (p == std::string::npos) return "./";
  else return filename.substr(0, p+1);
}

std::string add_suffix(const std::string &filename, const std::string &suffix) {
  auto stem = get_stem(filename);
  auto ext = get_ext(filename);
  return StringPrintf("%s_%s%s", stem.c_str(), suffix.c_str(), ext.c_str());
}

std::string add_index(const std::string &filename, int index,
                      const std::string &pattern) {
  auto stem = get_stem(filename);
  auto ext = get_ext(filename);
  return StringPrintf(StringPrintf("%s_%s%s", stem.c_str(),
                                   pattern.c_str(), ext.c_str()).c_str(),
                      index);
}

bool ListDir(const std::string &dirname, const std::string &ext,
             std::vector<std::string> *filenames) {
  DIR *pDIR;
	struct dirent *entry;
  pDIR = opendir(dirname.c_str());
  if(pDIR){
		while((entry = readdir(pDIR))){
			if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 &&
          (ext == "" || GetExt(entry->d_name) == ext))
      {
        filenames->push_back(entry->d_name);
      }
		}
    closedir(pDIR);
    return true;
	} else {
    return false;
  }
}

} // furry
