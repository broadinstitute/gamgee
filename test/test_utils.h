#ifndef gamgee_test_utils__guard
#define gamgee_test_utils__guard

#include <tuple>

/**
 * @brief test code for copy construction and copy assignment for any copy enabled object
 *
 * This tests copy construction, copy assignment and checks for self copy
 * assignment. The starting object is copied to three different objects using
 * all the above procedures. It returns a tuple with all three copied objects
 * for you to verify that the copy occurred and that it didn't affect the
 * original object. This is extremely useful when writing tests for a new class.
 * 
 * @tparam T any class that is copy constructible and copy assignable
 * @param original a simple object of class T to make the copies from
 * @return a tuple with all three copied objects
 */
template <class T>
std::tuple<T, T, T> check_copy_constructor (T& original) {
  auto obj2 = original; // copy construction
  auto obj3 = original; // copy construction
  obj2 = obj2;          // check self copy-assignment
  obj3 = obj2;          // copy assignment
  return std::make_tuple(original, obj2, obj3);
}

 /**
 * @brief test code for move construction and move assignment for any move enabled object
 *
 * This tests move construction, move assignment and checks for self move
 * assignment. The starting object is moved through three different objects using 
 * all the above procedures. It returns the final moved to object
 * for you to verify that the move occurred and that it still matches the 
 * original object. This is extremely useful when writing tests for a new class.
 *
 * @warning the original object passed in will be forcefully moved from, therefore in unusable state.
 * 
 * @tparam T any class that is move constructible and move assignable
 * @param original a simple object of class T to move from (it will be destroyed)
 * @return the last moved to object
 */
template <class T>
T check_move_constructor (T& original) {
  auto obj2 = std::move(original); // move construction
  auto obj3 = std::move(obj2);     // create a new object (we can't rely on default construction)
  obj2 = std::move(obj3);          // check move-assignment
  return obj2;
}



#endif
