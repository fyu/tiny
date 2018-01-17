#ifndef FURRY_COMMON_TEST_H_
#define FURRY_COMMON_TEST_H_

#include <functional>
#include <iostream>
#include <sstream>

#ifndef FURRY_CHECK
#define FURRY_CHECK 0
#endif

namespace furry {

#define FURRY_ASSERT(exp) \
  do { if (!(exp)) cerr << "Error: " #exp " is false on Line " << __LINE__ << " in File \"" << __FILE__ << "\"\n"; \
  } while (0)

#define F_ASSERT(exp) FURRY_ASSERT(exp)

namespace test {

class Message {
 public:
  Message() {}
  Message(const Message &msg) :
      ss_(msg.ss_.str()) {
  }
  template <typename T>
  Message& operator << (const T& value) {
    ss_ << value;
    return *this;
  }

  std::string GetString() const {
    return ss_.str();
  }
 private:
  std::stringstream ss_;
};

class TestResult {
 public:
  TestResult(bool result, int line, const std::string &file):
      result_(result), line_(line), file_(file) {}
  const TestResult& operator = (const Message &msg) const {
    if (!result_) {
      auto error_stream = &std::cout;
      *error_stream << "Error: " << file_ << " line: " << line_;
      auto str = msg.GetString();
      if (str != "")
        *error_stream << " assert: " << str;
      *error_stream << std::endl;
    }
    return *this;
  }
 private:
  bool result_;
  int line_;
  std::string file_;
};

} // test
} // furry

#define FCHECK(exp) \
  furry::check(exp, __LINE__, __FILE__) = furry::test::Message()

namespace furry {

inline test::TestResult check(std::function<bool ()> func,
                              int line, const std::string &file) {
#if FURRY_CHECK
  if (!func())
    return test::TestResult(false, line, file);
#endif
  return test::TestResult(true, line, file);
}

//void check(bool, const string &msg);
} // furry

#endif
