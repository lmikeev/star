/*
 *  html_helper.hpp
 *
 *  Created by Linar Mikeev (mikeev@cs.uni-saarland.de).
 *  Copyright (C) 2015 Saarland University. All rights reserved.
 *
 *
 *  This file is part of star.
 *
 *  star is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  star is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with star.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SOLVER_UTIL_HTML_HELPER_HPP_
#define SOLVER_UTIL_HTML_HELPER_HPP_

#include <iostream>

namespace solver {
namespace util {

class htmlHelper {
 public:
  static void write_beginning(std::ostream& os, char const* const title,
                              char const* const path) {
    os << "<html>" << std::endl;
    os << "<head>" << std::endl;
    os << "<title>" << title << "</title>" << std::endl;
    os << "<link rel='stylesheet' type='text/css' href='" << path
       << "/report.css'/>" << std::endl;
    os << "</head>" << std::endl;
    os << "<body>" << std::endl;
  }

  static void write_ending(std::ostream& os) {
    os << "</body>" << std::endl;
    os << "</html>" << std::endl;
  }

  static void write_begin_table(std::ostream& os) {
    os << "<table class='view-table'>" << std::endl;
  }

  static void write_end_table(std::ostream& os) {
    os << "</table>" << std::endl;
  }

  static void write_link(std::ostream& os, char const* const url,
                         char const* const text) {
    os << "<a href='" << url << "'>" << text << "</a>" << std::endl;
  }

  static void write_image(std::ostream& os, char const* const url,
                          char const* const text) {
    os << "<a href='" << url << "'><img class='plot' src='" << url << "' alt='"
       << text << "'></a>" << std::endl;
  }

  static void write_link(std::ostream& os, wchar_t const* const url,
                         char const* const text) {
    os << "<a href='" << url << "'>" << text << "</a>" << std::endl;
  }

  static void write_image(std::ostream& os, wchar_t const* const url,
                          char const* const text) {
    os << "<a href='" << url << "'><img class='plot' src='" << url << "' alt='"
       << text << "'></a>" << std::endl;
  }
};
}
}

#endif
