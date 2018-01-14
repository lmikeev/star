

#ifndef _GNUPLOT_PIPES_H_
#define _GNUPLOT_PIPES_H_

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <cstdlib>
#include <list>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)

#include <io.h>
#define GP_MAX_TMP_FILES 27
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)

#include <unistd.h>
#define GP_MAX_TMP_FILES 64
#else
#error unsupported or unknown operating system
#endif

class GnuplotException : public std::runtime_error {
 public:
  GnuplotException(const std::string &msg) : std::runtime_error(msg) {}
};

class Gnuplot {
 private:
  FILE *gnucmd;

  bool valid;

  bool two_dim;

  int nplots;

  std::string pstyle;

  std::string smooth;

  std::vector<std::string> tmpfile_list;

  static int tmpfile_num;

  static std::string m_sGNUPlotFileName;

  static std::string m_sGNUPlotPath;

  static std::string terminal_std;

  void init();

  std::string create_tmpfile(std::ofstream &tmp);

  static bool get_program_path();

  bool file_available(const std::string &filename);

  static bool file_exists(const std::string &filename, int mode = 0);

 public:
  static bool set_GNUPlotPath(const std::string &path);

  static void set_terminal_std(const std::string &type);

  Gnuplot(const std::string &style = "points");

  Gnuplot(const std::vector<double> &x, const std::string &title = "",
          const std::string &style = "points", const std::string &labelx = "x",
          const std::string &labely = "y");

  Gnuplot(const std::vector<double> &x, const std::vector<double> &y,
          const std::string &title = "", const std::string &style = "points",
          const std::string &labelx = "x", const std::string &labely = "y");

  Gnuplot(const std::vector<double> &x, const std::vector<double> &y,
          const std::vector<double> &z, const std::string &title = "",
          const std::string &style = "points", const std::string &labelx = "x",
          const std::string &labely = "y", const std::string &labelz = "z");

  ~Gnuplot();

  Gnuplot &cmd(const std::string &cmdstr);

  inline Gnuplot &operator<<(const std::string &cmdstr) {
    cmd(cmdstr);
    return (*this);
  }

  Gnuplot &showonscreen();

  Gnuplot &savetofigure(const std::string filename,
                        const std::string terminal = "ps");

  Gnuplot &set_style(const std::string &stylestr = "points");

  Gnuplot &set_smooth(const std::string &stylestr = "csplines");

  inline Gnuplot &unset_smooth() {
    smooth = "";
    return *this;
  };

  Gnuplot &set_pointsize(const double pointsize = 1.0);

  inline Gnuplot &set_grid() {
    cmd("set grid");
    return *this;
  };

  inline Gnuplot &unset_grid() {
    cmd("unset grid");
    return *this;
  };

  inline Gnuplot &set_multiplot() {
    cmd("set multiplot");
    return *this;
  };

  inline Gnuplot &unset_multiplot() {
    cmd("unset multiplot");
    return *this;
  };

  Gnuplot &set_samples(const int samples = 100);

  Gnuplot &set_isosamples(const int isolines = 10);

  Gnuplot &set_hidden3d() {
    cmd("set hidden3d");
    return *this;
  };

  inline Gnuplot &unset_hidden3d() {
    cmd("unset hidden3d");
    return *this;
  };

  Gnuplot &set_contour(const std::string &position = "base");

  inline Gnuplot &unset_contour() {
    cmd("unset contour");
    return *this;
  };

  inline Gnuplot &set_surface() {
    cmd("set surface");
    return *this;
  };

  inline Gnuplot &unset_surface() {
    cmd("unset surface");
    return *this;
  }

  Gnuplot &set_legend(const std::string &position = "default");

  inline Gnuplot &unset_legend() {
    cmd("unset key");
    return *this;
  }

  inline Gnuplot &set_title(const std::string &title = "") {
    std::string cmdstr;
    cmdstr = "set title \"";
    cmdstr += title;
    cmdstr += "\"";
    *this << cmdstr;
    return *this;
  }

  inline Gnuplot &unset_title() {
    this->set_title();
    return *this;
  }

  Gnuplot &set_ylabel(const std::string &label = "x");

  Gnuplot &set_xlabel(const std::string &label = "y");

  Gnuplot &set_zlabel(const std::string &label = "z");

  Gnuplot &set_xrange(const double iFrom, const double iTo);

  Gnuplot &set_yrange(const double iFrom, const double iTo);

  Gnuplot &set_zrange(const double iFrom, const double iTo);

  inline Gnuplot &set_xautoscale() {
    cmd("set xrange restore");
    cmd("set autoscale x");
    return *this;
  };

  inline Gnuplot &set_yautoscale() {
    cmd("set yrange restore");
    cmd("set autoscale y");
    return *this;
  };

  inline Gnuplot &set_zautoscale() {
    cmd("set zrange restore");
    cmd("set autoscale z");
    return *this;
  };

  Gnuplot &set_xlogscale(const double base = 10);

  Gnuplot &set_ylogscale(const double base = 10);

  Gnuplot &set_zlogscale(const double base = 10);

  inline Gnuplot &unset_xlogscale() {
    cmd("unset logscale x");
    return *this;
  };

  inline Gnuplot &unset_ylogscale() {
    cmd("unset logscale y");
    return *this;
  };

  inline Gnuplot &unset_zlogscale() {
    cmd("unset logscale z");
    return *this;
  };

  Gnuplot &set_cbrange(const double iFrom, const double iTo);

  Gnuplot &plotfile_x(const std::string &filename,
                      const unsigned int column = 1,
                      const std::string &title = "");

  template <typename X>
  Gnuplot &plot_x(const X &x, const std::string &title = "");

  Gnuplot &plotfile_xy(const std::string &filename,
                       const unsigned int column_x = 1,
                       const unsigned int column_y = 2,
                       const std::string &title = "");

  template <typename X, typename Y>
  Gnuplot &plot_xy(const X &x, const Y &y, const std::string &title = "");

  Gnuplot &plotfile_xy_err(const std::string &filename,
                           const unsigned int column_x = 1,
                           const unsigned int column_y = 2,
                           const unsigned int column_dy = 3,
                           const std::string &title = "");

  template <typename X, typename Y, typename E>
  Gnuplot &plot_xy_err(const X &x, const Y &y, const E &dy,
                       const std::string &title = "");

  Gnuplot &plotfile_xyz(const std::string &filename,
                        const unsigned int column_x = 1,
                        const unsigned int column_y = 2,
                        const unsigned int column_z = 3,
                        const std::string &title = "");

  template <typename X, typename Y, typename Z>
  Gnuplot &plot_xyz(const X &x, const Y &y, const Z &z,
                    const std::string &title = "");

  Gnuplot &plot_slope(const double a, const double b,
                      const std::string &title = "");

  Gnuplot &plot_equation(const std::string &equation,
                         const std::string &title = "");

  Gnuplot &plot_equation3d(const std::string &equation,
                           const std::string &title = "");

  Gnuplot &plot_image(const unsigned char *ucPicBuf, const unsigned int iWidth,
                      const unsigned int iHeight,
                      const std::string &title = "");

  inline Gnuplot &replot(void) {
    if (nplots > 0) cmd("replot");
    return *this;
  };

  Gnuplot &reset_plot();

  Gnuplot &reset_all();

  void remove_tmpfiles();

  inline bool is_valid() { return (valid); };
};

int Gnuplot::tmpfile_num = 0;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
std::string Gnuplot::m_sGNUPlotFileName = "pgnuplot.exe";
std::string Gnuplot::m_sGNUPlotPath = "C:/program files/gnuplot/bin/";
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
std::string Gnuplot::m_sGNUPlotFileName = "gnuplot";
std::string Gnuplot::m_sGNUPlotPath = "/usr/local/bin/";
#endif

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
std::string Gnuplot::terminal_std = "windows";
#elif(defined(unix) || defined(__unix) || defined(__unix__)) && \
    !defined(__APPLE__)
std::string Gnuplot::terminal_std = "x11";
#elif defined(__APPLE__)
std::string Gnuplot::terminal_std = "aqua";
#endif

inline Gnuplot::Gnuplot(const std::string &style)
    : gnucmd(NULL),
      valid(false),
      two_dim(false),
      nplots(0)

{
  init();
  set_style(style);
}

inline Gnuplot::Gnuplot(const std::vector<double> &x, const std::string &title,
                        const std::string &style, const std::string &labelx,
                        const std::string &labely)
    : gnucmd(NULL), valid(false), two_dim(false), nplots(0) {
  init();

  set_style(style);
  set_xlabel(labelx);
  set_ylabel(labely);

  plot_x(x, title);
}

inline Gnuplot::Gnuplot(const std::vector<double> &x,
                        const std::vector<double> &y, const std::string &title,
                        const std::string &style, const std::string &labelx,
                        const std::string &labely)
    : gnucmd(NULL), valid(false), two_dim(false), nplots(0) {
  init();

  set_style(style);
  set_xlabel(labelx);
  set_ylabel(labely);

  plot_xy(x, y, title);
}

inline Gnuplot::Gnuplot(const std::vector<double> &x,
                        const std::vector<double> &y,
                        const std::vector<double> &z, const std::string &title,
                        const std::string &style, const std::string &labelx,
                        const std::string &labely, const std::string &labelz)
    : gnucmd(NULL), valid(false), two_dim(false), nplots(0) {
  init();

  set_style(style);
  set_xlabel(labelx);
  set_ylabel(labely);
  set_zlabel(labelz);

  plot_xyz(x, y, z, title);
}

template <typename X>
Gnuplot &Gnuplot::plot_x(const X &x, const std::string &title) {
  if (x.size() == 0) {
    throw GnuplotException("std::vector too small");
    return *this;
  }

  std::ofstream tmp;
  std::string name = create_tmpfile(tmp);
  if (name == "") return *this;

  for (unsigned int i = 0; i < x.size(); i++) tmp << x[i] << std::endl;

  tmp.flush();
  tmp.close();

  plotfile_x(name, 1, title);

  return *this;
}

template <typename X, typename Y>
Gnuplot &Gnuplot::plot_xy(const X &x, const Y &y, const std::string &title) {
  if (x.size() == 0 || y.size() == 0) {
    throw GnuplotException("std::vectors too small");
    return *this;
  }

  if (x.size() != y.size()) {
    throw GnuplotException("Length of the std::vectors differs");
    return *this;
  }

  std::ofstream tmp;
  std::string name = create_tmpfile(tmp);
  if (name == "") return *this;

  for (unsigned int i = 0; i < x.size(); i++)
    tmp << x[i] << " " << y[i] << std::endl;

  tmp.flush();
  tmp.close();

  plotfile_xy(name, 1, 2, title);

  return *this;
}

template <typename X, typename Y, typename E>
Gnuplot &Gnuplot::plot_xy_err(const X &x, const Y &y, const E &dy,
                              const std::string &title) {
  if (x.size() == 0 || y.size() == 0 || dy.size() == 0) {
    throw GnuplotException("std::vectors too small");
    return *this;
  }

  if (x.size() != y.size() || y.size() != dy.size()) {
    throw GnuplotException("Length of the std::vectors differs");
    return *this;
  }

  std::ofstream tmp;
  std::string name = create_tmpfile(tmp);
  if (name == "") return *this;

  for (unsigned int i = 0; i < x.size(); i++)
    tmp << x[i] << " " << y[i] << " " << dy[i] << std::endl;

  tmp.flush();
  tmp.close();

  plotfile_xy_err(name, 1, 2, 3, title);

  return *this;
}

template <typename X, typename Y, typename Z>
Gnuplot &Gnuplot::plot_xyz(const X &x, const Y &y, const Z &z,
                           const std::string &title) {
  if (x.size() == 0 || y.size() == 0 || z.size() == 0) {
    throw GnuplotException("std::vectors too small");
    return *this;
  }

  if (x.size() != y.size() || x.size() != z.size()) {
    throw GnuplotException("Length of the std::vectors differs");
    return *this;
  }

  std::ofstream tmp;
  std::string name = create_tmpfile(tmp);
  if (name == "") return *this;

  for (unsigned int i = 0; i < x.size(); i++)
    tmp << x[i] << " " << y[i] << " " << z[i] << std::endl;

  tmp.flush();
  tmp.close();

  plotfile_xyz(name, 1, 2, 3, title);

  return *this;
}

bool Gnuplot::set_GNUPlotPath(const std::string &path) {
  std::string tmp = path + "/" + Gnuplot::m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  if (Gnuplot::file_exists(tmp, 0))
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  if (Gnuplot::file_exists(tmp, 1))
#endif
  {
    Gnuplot::m_sGNUPlotPath = path;
    return true;
  } else {
    Gnuplot::m_sGNUPlotPath.clear();
    return false;
  }
}

void Gnuplot::set_terminal_std(const std::string &type) {
#if defined(unix) || defined(__unix) || defined(__unix__) || defined(__APPLE__)
  if (type.find("x11") != std::string::npos && getenv("DISPLAY") == NULL) {
    throw GnuplotException("Can't find DISPLAY variable");
  }
#endif

  Gnuplot::terminal_std = type;
  return;
}

template <typename Container>
void stringtok(Container &container, std::string const &in,
               const char *const delimiters = " \t\n") {
  const std::string::size_type len = in.length();
  std::string::size_type i = 0;

  while (i < len) {
    i = in.find_first_not_of(delimiters, i);

    if (i == std::string::npos) return;

    std::string::size_type j = in.find_first_of(delimiters, i);

    if (j == std::string::npos) {
      container.push_back(in.substr(i));
      return;
    } else
      container.push_back(in.substr(i, j - i));

    i = j + 1;
  }

  return;
}

Gnuplot::~Gnuplot() {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  if (_pclose(gnucmd) == -1)
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  if (pclose(gnucmd) == -1)
#endif
    throw GnuplotException("Problem closing communication to gnuplot");
}

Gnuplot &Gnuplot::reset_plot() {
  nplots = 0;

  return *this;
}

Gnuplot &Gnuplot::reset_all() {
  nplots = 0;
  cmd("reset");
  cmd("clear");
  pstyle = "points";
  smooth = "";
  showonscreen();

  return *this;
}

Gnuplot &Gnuplot::set_style(const std::string &stylestr) {
  if (stylestr.find("lines") == std::string::npos &&
      stylestr.find("points") == std::string::npos &&
      stylestr.find("linespoints") == std::string::npos &&
      stylestr.find("impulses") == std::string::npos &&
      stylestr.find("dots") == std::string::npos &&
      stylestr.find("steps") == std::string::npos &&
      stylestr.find("fsteps") == std::string::npos &&
      stylestr.find("histeps") == std::string::npos &&
      stylestr.find("boxes") == std::string::npos &&
      stylestr.find("filledcurves") == std::string::npos &&
      stylestr.find("histograms") == std::string::npos)

  {
    pstyle = std::string("points");
  } else {
    pstyle = stylestr;
  }

  return *this;
}

Gnuplot &Gnuplot::set_smooth(const std::string &stylestr) {
  if (stylestr.find("unique") == std::string::npos &&
      stylestr.find("frequency") == std::string::npos &&
      stylestr.find("csplines") == std::string::npos &&
      stylestr.find("acsplines") == std::string::npos &&
      stylestr.find("bezier") == std::string::npos &&
      stylestr.find("sbezier") == std::string::npos) {
    smooth = "";
  } else {
    smooth = stylestr;
  }

  return *this;
}

Gnuplot &Gnuplot::showonscreen() {
  cmd("set output");
  cmd("set terminal " + Gnuplot::terminal_std);

  return *this;
}

Gnuplot &Gnuplot::savetofigure(const std::string filename,
                               const std::string terminal) {
  std::ostringstream cmdstr;
  cmdstr << "set terminal " << terminal;
  cmd(cmdstr.str());

  cmdstr.str("");
  cmdstr << "set output \"" << filename << "\"";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_legend(const std::string &position) {
  std::ostringstream cmdstr;
  cmdstr << "set key " << position;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_xlogscale(const double base) {
  std::ostringstream cmdstr;

  cmdstr << "set logscale x " << base;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_ylogscale(const double base) {
  std::ostringstream cmdstr;

  cmdstr << "set logscale y " << base;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_zlogscale(const double base) {
  std::ostringstream cmdstr;

  cmdstr << "set logscale z " << base;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_pointsize(const double pointsize) {
  std::ostringstream cmdstr;
  cmdstr << "set pointsize " << pointsize;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_samples(const int samples) {
  std::ostringstream cmdstr;
  cmdstr << "set samples " << samples;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_isosamples(const int isolines) {
  std::ostringstream cmdstr;
  cmdstr << "set isosamples " << isolines;
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_contour(const std::string &position) {
  if (position.find("base") == std::string::npos &&
      position.find("surface") == std::string::npos &&
      position.find("both") == std::string::npos) {
    cmd("set contour base");
  } else {
    cmd("set contour " + position);
  }

  return *this;
}

Gnuplot &Gnuplot::set_xlabel(const std::string &label) {
  std::ostringstream cmdstr;

  cmdstr << "set xlabel \"" << label << "\"";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_ylabel(const std::string &label) {
  std::ostringstream cmdstr;

  cmdstr << "set ylabel \"" << label << "\"";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_zlabel(const std::string &label) {
  std::ostringstream cmdstr;

  cmdstr << "set zlabel \"" << label << "\"";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_xrange(const double iFrom, const double iTo) {
  std::ostringstream cmdstr;

  cmdstr << "set xrange[" << iFrom << ":" << iTo << "]";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_yrange(const double iFrom, const double iTo) {
  std::ostringstream cmdstr;

  cmdstr << "set yrange[" << iFrom << ":" << iTo << "]";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_zrange(const double iFrom, const double iTo) {
  std::ostringstream cmdstr;

  cmdstr << "set zrange[" << iFrom << ":" << iTo << "]";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::set_cbrange(const double iFrom, const double iTo) {
  std::ostringstream cmdstr;

  cmdstr << "set cbrange[" << iFrom << ":" << iTo << "]";
  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plot_slope(const double a, const double b,
                             const std::string &title) {
  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  cmdstr << a << " * x + " << b << " title \"";

  if (title == "")
    cmdstr << "f(x) = " << a << " * x + " << b;
  else
    cmdstr << title;

  cmdstr << "\" with " << pstyle;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plot_equation(const std::string &equation,
                                const std::string &title) {
  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  cmdstr << equation;

  if (title == "")
    cmdstr << " notitle";
  else
    cmdstr << " title \"" << title << "\"";

  cmdstr << " with " << pstyle;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plot_equation3d(const std::string &equation,
                                  const std::string &title) {
  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == false)
    cmdstr << "replot ";
  else
    cmdstr << "splot ";

  cmdstr << equation << " title \"";

  if (title == "")
    cmdstr << "f(x,y) = " << equation;
  else
    cmdstr << title;

  cmdstr << "\" with " << pstyle;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plotfile_x(const std::string &filename,
                             const unsigned int column,
                             const std::string &title) {
  file_available(filename);

  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  cmdstr << "\"" << filename << "\" using " << column;

  if (title == "")
    cmdstr << " notitle ";
  else
    cmdstr << " title \"" << title << "\" ";

  if (smooth == "")
    cmdstr << "with " << pstyle;
  else
    cmdstr << "smooth " << smooth;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plotfile_xy(const std::string &filename,
                              const unsigned int column_x,
                              const unsigned int column_y,
                              const std::string &title) {
  file_available(filename);

  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  cmdstr << "\"" << filename << "\" using " << column_x << ":" << column_y;

  if (title == "")
    cmdstr << " notitle ";
  else
    cmdstr << " title \"" << title << "\" ";

  if (smooth == "")
    cmdstr << "with " << pstyle;
  else
    cmdstr << "smooth " << smooth;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plotfile_xy_err(const std::string &filename,
                                  const unsigned int column_x,
                                  const unsigned int column_y,
                                  const unsigned int column_dy,
                                  const std::string &title) {
  file_available(filename);

  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  cmdstr << "\"" << filename << "\" using " << column_x << ":" << column_y
         << ":" << column_dy << " with errorbars ";

  if (title == "")
    cmdstr << " notitle ";
  else
    cmdstr << " title \"" << title << "\" ";

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plotfile_xyz(const std::string &filename,
                               const unsigned int column_x,
                               const unsigned int column_y,
                               const unsigned int column_z,
                               const std::string &title) {
  file_available(filename);

  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == false)
    cmdstr << "replot ";
  else
    cmdstr << "splot ";

  cmdstr << "\"" << filename << "\" using " << column_x << ":" << column_y
         << ":" << column_z;

  if (title == "")
    cmdstr << " notitle with " << pstyle;
  else
    cmdstr << " title \"" << title << "\" with " << pstyle;

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::plot_image(const unsigned char *ucPicBuf,
                             const unsigned int iWidth,
                             const unsigned int iHeight,
                             const std::string &title) {
  std::ofstream tmp;
  std::string name = create_tmpfile(tmp);
  if (name == "") return *this;

  int iIndex = 0;
  for (unsigned int iRow = 0; iRow < iHeight; iRow++) {
    for (unsigned int iColumn = 0; iColumn < iWidth; iColumn++) {
      tmp << iColumn << " " << iRow << " "
          << static_cast<float>(ucPicBuf[iIndex++]) << std::endl;
    }
  }

  tmp.flush();
  tmp.close();

  std::ostringstream cmdstr;

  if (nplots > 0 && two_dim == true)
    cmdstr << "replot ";
  else
    cmdstr << "plot ";

  if (title == "")
    cmdstr << "\"" << name << "\" with image";
  else
    cmdstr << "\"" << name << "\" title \"" << title << "\" with image";

  cmd(cmdstr.str());

  return *this;
}

Gnuplot &Gnuplot::cmd(const std::string &cmdstr) {
  if (!(valid)) {
    return *this;
  }

  fputs((cmdstr + "\n").c_str(), gnucmd);

  fflush(gnucmd);

  if (cmdstr.find("replot") != std::string::npos) {
    return *this;
  } else if (cmdstr.find("splot") != std::string::npos) {
    two_dim = false;
    nplots++;
  } else if (cmdstr.find("plot") != std::string::npos) {
    two_dim = true;
    nplots++;
  }

  return *this;
}

void Gnuplot::init() {
#if (defined(unix) || defined(__unix) || defined(__unix__)) && \
    !defined(__APPLE__)
  if (getenv("DISPLAY") == NULL) {
    valid = false;
    throw GnuplotException("Can't find DISPLAY variable");
  }
#endif

  if (!Gnuplot::get_program_path()) {
    valid = false;
    throw GnuplotException("Can't find gnuplot");
  }

  std::string tmp = Gnuplot::m_sGNUPlotPath + "/" + Gnuplot::m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  gnucmd = _popen(tmp.c_str(), "w");
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  gnucmd = popen(tmp.c_str(), "w");
#endif

  if (!gnucmd) {
    valid = false;
    throw GnuplotException("Couldn't open connection to gnuplot");
  }

  nplots = 0;
  valid = true;
  smooth = "";

  showonscreen();

  return;
}

bool Gnuplot::get_program_path() {
  std::string tmp = Gnuplot::m_sGNUPlotPath + "/" + Gnuplot::m_sGNUPlotFileName;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  if (Gnuplot::file_exists(tmp, 0))
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  if (Gnuplot::file_exists(tmp, 1))
#endif
  {
    return true;
  }

  char *path;

  path = getenv("PATH");

  if (path == NULL) {
    throw GnuplotException("Path is not set");
    return false;
  } else {
    std::list<std::string> ls;

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
    stringtok(ls, path, ";");
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
    stringtok(ls, path, ":");
#endif

    for (std::list<std::string>::const_iterator i = ls.begin(); i != ls.end();
         ++i) {
      tmp = (*i) + "/" + Gnuplot::m_sGNUPlotFileName;
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
      if (Gnuplot::file_exists(tmp, 0))
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
      if (Gnuplot::file_exists(tmp, 1))
#endif
      {
        Gnuplot::m_sGNUPlotPath = *i;
        return true;
      }
    }

    tmp = "Can't find gnuplot neither in PATH nor in \"" +
          Gnuplot::m_sGNUPlotPath + "\"";
    throw GnuplotException(tmp);

    Gnuplot::m_sGNUPlotPath = "";
    return false;
  }
}

bool Gnuplot::file_exists(const std::string &filename, int mode) {
  if (mode < 0 || mode > 7) {
    throw std::runtime_error(
        "In function \"Gnuplot::file_exists\": mode\
                has to be an integer between 0 and 7");
    return false;
  }

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  if (_access(filename.c_str(), mode) == 0)
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  if (access(filename.c_str(), mode) == 0)
#endif
  {
    return true;
  } else {
    return false;
  }
}

bool Gnuplot::file_available(const std::string &filename) {
  std::ostringstream except;
  if (Gnuplot::file_exists(filename, 0)) {
    if (!(Gnuplot::file_exists(filename, 4))) {
      except << "No read permission for File \"" << filename << "\"";
      throw GnuplotException(except.str());
      return false;
    }
  } else {
    except << "File \"" << filename << "\" does not exist";
    throw GnuplotException(except.str());
    return false;
  }
  return false;
}

std::string Gnuplot::create_tmpfile(std::ofstream &tmp) {
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  char name[] = "gnuplotiXXXXXX";
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  char name[] = "/tmp/gnuplotiXXXXXX";
#endif

  if (Gnuplot::tmpfile_num == GP_MAX_TMP_FILES - 1) {
    std::ostringstream except;
    except << "Maximum number of temporary files reached (" << GP_MAX_TMP_FILES
           << "): cannot open more files" << std::endl;

    throw GnuplotException(except.str());
    return "";
  }

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || \
    defined(__TOS_WIN__)
  if (_mktemp(name) == NULL)
#elif defined(unix) || defined(__unix) || defined(__unix__) || \
    defined(__APPLE__)
  if (mkstemp(name) == -1)
#endif
  {
    std::ostringstream except;
    except << "Cannot create temporary file \"" << name << "\"";
    throw GnuplotException(except.str());
  }

  tmp.open(name);
  if (tmp.bad()) {
    std::ostringstream except;
    except << "Cannot create temporary file \"" << name << "\"";
    throw GnuplotException(except.str());
  }

  tmpfile_list.push_back(name);
  Gnuplot::tmpfile_num++;

  return name;
}

void Gnuplot::remove_tmpfiles() {
  if ((tmpfile_list).size() > 0) {
    for (unsigned int i = 0; i < tmpfile_list.size(); i++) {
      if (remove(tmpfile_list[i].c_str()) != 0) {
        std::ostringstream except;
        except << "Cannot remove temporary file \"" << tmpfile_list[i] << "\"";
        throw GnuplotException(except.str());
      }
    }

    Gnuplot::tmpfile_num -= tmpfile_list.size();
  }
}
#endif
