#ifndef _FURRY_CORE_FMETA_H_
#define _FURRY_CORE_FMETA_H_

namespace furry
{

#define APPLY(f, x) f::apply<x>::type

template <class f, class x>
struct apply : f::template apply<x>
{
};

template <class f, class x1, class x2>
struct apply2 : f::template apply<x1>::type::template apply<x2>
{
};

struct first_
{
  template <class x, class ... xs>
  struct f
  {
    typedef x type;
  };
};

struct second_
{
  template <class x, class ... xs>
  struct f : first_::template f<xs...> {};
};

/*struct true_
{
  template <typename x1>
  struct f
  {
    typedef struct
    {
      template <typename x2>
      struct f
      {
        typedef x1 type;
      };
    } type;
  };
};

struct false_
{
  template <typename x1>
  struct f
  {
    typedef struct
    {
      template <typename x2>
      struct f
      {
        typedef x2 type;
      };
    } type;
  };
};

struct if
{
  template <class b>
  struct f
  {
    typedef struct
    {
      template <class x1>
      struct f
      {
        typedef struct
        {
          template <class x2>
          struct f : f2<b
        }
      }
    }
  }
}*/

struct true_
{
  template <class x1, class x2>
  struct f
  {
    typedef x1 type;
  };
};

struct false_
{
  template <class x1, class x2>
  struct f
  {
    typedef x2 type;
  };
};

template <bool b>
struct bool_ : true_ {};

template <>
struct bool_<false> : false_ {};

struct if_
{
  template <class b, class x1, class x2>
  struct f
  {
    typedef typename b::template f<x1, x2>::type type;
  };
};

struct and_
{
  template <class b1, class b2>
  struct f : b1::template f<b2, false_> {};
};

struct or_
{
  template <class b1, class b2>
  struct f : b1::template f<true_, b2> {};
};

struct not_
{
  template <class b>
  struct f : b::template f<false_, true_> {};
};

} // furry

#endif // _FURRY_CORE_FMETA_H_
