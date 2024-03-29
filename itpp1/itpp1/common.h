#ifndef HAVE_COMMON_H
#define HAVE_COMMON_H

// #define USE_OPENCL 0

// #define HAVE_RTLSDR 0

// #define HAVE_HACKRF 0

// #define HAVE_BLADERF 0

// This is filled in by cmake
// #define MAJOR_VERSION @CellSearch_MAJOR_VERSION@
// #define MINOR_VERSION @CellSearch_MINOR_VERSION@
// #define PATCH_LEVEL @CellSearch_PATCH_LEVEL@
// #define BUILD_TYPE "@CMAKE_BUILD_TYPE@"

// Typedefs
typedef char int8;
typedef unsigned char uint8;
typedef short int16;
typedef unsigned short uint16;
typedef int int32;
typedef unsigned int uint32;
typedef long int int64;
typedef unsigned long int uint64;
typedef unsigned int uint;
// Prevent int8 and uint8 from printing out as characters instead of integers.
inline std::ostream & operator<< (
  std::ostream & os,
  const uint8 & c
) {
  os << ((unsigned int)c);
  return os;
}
inline std::ostream & operator<< (
  std::ostream & os,
  const int8 & c
) {
  os << ((int)c);
  return os;
}

// 0 indicates minimum amount of status messages to cout
// 1 indicates normal amount of status messages to cout
// 2 indicates maximum amount of status messages to cout
extern uint8 verbosity;

// complex<float> 2d/3d vectors
typedef std::vector < std::vector < std::vector < std::complex < float > > > > vcf3d;
typedef std::vector < std::vector < std::complex < float > > > vcf2d;
typedef std::vector < std::vector < std::vector < float > > > vf3d;
typedef std::vector < std::vector < float > > vf2d;

// Some enums must be enclosed in their own namespace because their
// names conflict with each other and also with ITPP declared enums.

namespace dev_type_t {
  enum dev_type_t { UNKNOWN = -4321, RTLSDR=9832, HACKRF=432134, BLADERF=94703 };
}

namespace cp_type_t {
  enum cp_type_t { UNKNOWN = 0, NORMAL, EXTENDED };
}
inline std::ostream & operator<< (
  std::ostream & os,
  const cp_type_t::cp_type_t & c
) {
  switch (c) {
    case cp_type_t::UNKNOWN: os << "UNKNOWN"; break;
    case cp_type_t::NORMAL: os << "NORMAL"; break;
    case cp_type_t::EXTENDED: os << "EXTENDED"; break;
    default: os << "???"; break;
  }
  return os;
}
enum crc_t { CRC8, CRC16, CRC24A, CRC24B };
namespace phich_duration_t {
  enum phich_duration_t {UNKNOWN = 0, NORMAL, EXTENDED};
}
inline std::ostream & operator<< (
  std::ostream & os,
  const phich_duration_t::phich_duration_t & c
) {
  switch (c) {
    case phich_duration_t::UNKNOWN: os << "UNKNOWN"; break;
    case phich_duration_t::NORMAL: os << "NORMAL"; break;
    case phich_duration_t::EXTENDED: os << "EXTENDED"; break;
    default: os << "???"; break;
  }
  return os;
}
namespace phich_resource_t {
  enum phich_resource_t {UNKNOWN = 0, oneSixth, half, one, two};
}
inline std::ostream & operator<< (
  std::ostream & os,
  const phich_resource_t::phich_resource_t & c
) {
  switch (c) {
    case phich_resource_t::UNKNOWN: os << "UNKNOWN"; break;
    case phich_resource_t::oneSixth: os << "oneSixth"; break;
    case phich_resource_t::half: os << "half"; break;
    case phich_resource_t::one: os << "one"; break;
    case phich_resource_t::two: os << "two"; break;
    default: os << "???"; break;
  }
  return os;
}
namespace modulation_t {
  enum modulation_t {QAM, QAM16, QAM64};
}

// Class to contain all the information detected about a cell.
class Cell {
  public:
    double fc_requested;
    double fc_programmed;
    double pss_pow;
    int32 ind;
    double freq;
    int8 n_id_2;
    double k_factor;

    int16 n_id_1;
    int8 duplex_mode;
    cp_type_t::cp_type_t cp_type;
    double frame_start;
    double freq_fine;

    double freq_superfine;

    int8 n_ports;
    int8 n_rb_dl;
    phich_duration_t::phich_duration_t phich_duration;
    phich_resource_t::phich_resource_t phich_resource;
    int16 sfn;
    // Member functions
    // Constructors
    Cell();
    // Misc
    int16 const n_id_cell() const;
    int8 const n_symb_dl() const;
  private:
};

// Allow for easy printing of a 'Cell'
std::ostream & operator<< (
  std::ostream & os,
  const Cell & c
);

// The delay is necessary to prevent a sefault happening on exit.
#define ABORT(A) endwin(); usleep(100); exit(A)

#endif

