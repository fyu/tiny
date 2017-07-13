#ifndef _FURRY_CORE_META_H_
#define _FURRY_CORE_MAT_H_

namespace furry
{

#define META_VALUE(value) typename value##::type

template <typename T1, typename T2>
struct true_
{
  typedef T1 type;
};

template <typename T1, typename T2>
struct false_
{
  typedef T2 type;
};

template <bool>
struct to_bool_
{
  template <typename T1, typename T2>
  using type = typename true_<T1, T2>::type;
};

template <>
struct to_bool_<false>
{
  typedef false_ type;
};

#define TO_BOOL(b) typename to_bool_<b>::type

//template <template <typename B1T1, typename B1T2> class B1,
//        template <typename B2T1, typename B2T2> class B2>

template <template <typename BT1, typename BT2> class B,
          typename T1,
          typename T2>
struct if_
{
  typedef typename B<T1, T2>::type type;
};

} // furry

#endif // _FURRY_CORE_MAT_H_
