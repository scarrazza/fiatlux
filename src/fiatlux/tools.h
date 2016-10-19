
// FiatLux 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch

#pragma once

#include <iostream>
#include <string>
using std::string;
using std::cout;
using std::endl;

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

}
