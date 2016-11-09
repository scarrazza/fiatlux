
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <iostream>
#include <string>
#include <sstream>
using std::string;
using std::cout;
using std::endl;
using std::stringstream;

namespace fiatlux {

  /**
   * @brief Just a simple print layout function.
   *
   * @param tag the function name.
   * @param text the message content.
   */
  inline void info(string const& tag, string const& text)
  {
    cout << "[" << tag << "] " << text << endl;
  }

  inline string infos(string const& tag, string const& text)
  {
    stringstream ss("");
    ss << "[" << tag << "] " << text << endl;
    return ss.str();
  }

  /**
   * @brief The runtime_exception class
   */
  class runtime_exception: public std::runtime_error
  {
  public:
    runtime_exception(string const& tag, string const& what):
      std::runtime_error(infos(tag, what)) {}
  };

  /**
   * @brief swap method for doubles
   */
  inline void swap(double* x, double* y)
  {
    double tmp;
    tmp = *y;
    *y = *x;
    *x = tmp;
  }  
}
