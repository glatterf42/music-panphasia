#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <string>

bool isFloat(std::string myString) {
  std::istringstream iss(myString);
  double f;
  // iss >> std::noskipws >> f; // noskipws considers leading whitespace invalid
  // Check the entire string was consumed and if either failbit or badbit is set
  iss >> f;
  return iss.eof() && !iss.fail();
}

template <typename real_t>
void read_points_from_file(std::string fname, float vfac_,
                           std::vector<real_t> &p) {

  int num_columns;

  std::ifstream ifs(fname.c_str());
  if (!ifs) {
    fprintf(stderr,
            "region_ellipsoid_plugin::read_points_from_file : Could not open "
            "file \'%s\'",
            fname.c_str());
    throw std::runtime_error("region_ellipsoid_plugin::read_points_from_file : "
                             "cannot open point file.");
  }

  int colcount = 0, colcount1 = 0, row = 0;
  p.clear();

  while (ifs) {
    std::string s;
    if (!getline(ifs, s))
      break;
    //std::cerr << s << std::endl;
    std::replace(s.begin(), s.end(), '\t', ' ');
    std::replace(s.begin(), s.end(), '\r', ' ');
    std::stringstream ss(s);
    colcount1 = 0;
    while (ss) {
      if (!getline(ss, s, ' '))
        break;
      if (!isFloat(s))
        continue;
      p.push_back(strtod(s.c_str(), NULL));

      if (row == 0)
        colcount++;
      else
        colcount1++;
    }
    ++row;

    if (row > 1 && colcount != colcount1)
      fprintf(stderr, "error on line %d of input file", row);

    // std::cout << std::endl;
  }

  printf("region point file appears to contain %d columns", colcount);

  if (p.size() % 3 != 0 && p.size() % 6 != 0) {
    for (auto x : p)
      std::cerr << x << "\n";
    fprintf(stderr,
            "Region point file \'%s\' does not contain triplets (%d elems)",
            fname.c_str(), p.size());
    throw std::runtime_error("region_ellipsoid_plugin::read_points_from_file : "
                             "file does not contain triplets.");
  }

  double x0[3] = {p[0], p[1], p[2]}, dx;

  if (colcount == 3) {
    // only positions are given

    for (size_t i = 3; i < p.size(); i += 3) {
      for (size_t j = 0; j < 3; ++j) {
        dx = p[i + j] - x0[j];
        if (dx < -0.5)
          dx += 1.0;
        else if (dx > 0.5)
          dx -= 1.0;
        p[i + j] = x0[j] + dx;
      }
    }
  } else if (colcount == 6) {
    // positions and velocities are given

    //... include the velocties to unapply Zeldovich approx.

    for (size_t j = 3; j < 6; ++j) {
      dx = (p[j - 3] - p[j] / vfac_) - x0[j - 3];
      if (dx < -0.5)
        dx += 1.0;
      else if (dx > 0.5)
        dx -= 1.0;
      p[j] = x0[j - 3] + dx;
    }

    for (size_t i = 6; i < p.size(); i += 6) {
      for (size_t j = 0; j < 3; ++j) {
        dx = p[i + j] - x0[j];
        if (dx < -0.5)
          dx += 1.0;
        else if (dx > 0.5)
          dx -= 1.0;
        p[i + j] = x0[j] + dx;
      }

      for (size_t j = 3; j < 6; ++j) {
        dx = (p[i + j - 3] - p[i + j] / vfac_) - x0[j - 3];
        if (dx < -0.5)
          dx += 1.0;
        else if (dx > 0.5)
          dx -= 1.0;
        p[i + j] = x0[j - 3] + dx;
      }
    }
  } else
    fprintf(stderr, "Problem interpreting the region point file \'%s\'",
            fname.c_str());

  num_columns = colcount;
  fprintf(stderr,"num_columns = %d\n\n",num_columns);
}



int main( int argc, char**argv )
{
  std::vector<float> pts;
  read_points_from_file<float>( argv[1], 1.0, pts );

  return 1;
}
