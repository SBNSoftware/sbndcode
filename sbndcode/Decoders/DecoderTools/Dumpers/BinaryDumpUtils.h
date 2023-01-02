/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/BinaryDumpUtils.h
 * @brief  Functions to dump the content of binary data chunks to console
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_BINARYDUMPUTILS_H
#define SBNDCODE_DECODERS_DECODERTOOLS_BINARYDUMPUTILS_H

// C/C++ standard libraries
#include <ostream>
#include <iomanip> // std::setw()
#include <utility> // std::move
#include <type_traits> // std::is_integral_v
#include <cstddef> // std::size_t, std::ptrdiff_t


// -----------------------------------------------------------------------------
namespace sbnd::ns::util::details {
  
  // --- BEGIN -- Wrapping objects ---------------------------------------------
  /**
   * @name Wrapping binary objects
   * 
   * The objects in this section are implementation details.
   * Interface functions convert objects into these helpers, and specific
   * output operators (`<<`) will manage their output to stream.
   * 
   * This is a mechanism similar to C++ STL I/O manipulators with arguments.
   */
  /// @{
  
  /**
   * @brief An object wrapping some data (copy), with a tag type.
   * @tparam Tag a type to make a template class instance unique
   * @tparam T type of the data wrapped
   * @tparam Bits size of data in `T` in bits
   * 
   * This template is used to implement other wrappers.
   * 
   * The role of `Tag` is just to allow distinguishing template instances with
   * the same `T` (see for example `BinObj` and `HexObj`).
   */
  template <typename Tag, typename T, unsigned int const Bits = (8*sizeof(T))>
  struct BitObjHolder {
    T const data;
    static constexpr unsigned int bits { Bits };
    constexpr BitObjHolder(T data): data(data) {}
  }; // BitObjHolder<>
  
  
  struct BinObjTag {}; ///< Tag object for `BinObj`.
  
  /// Holder for data to be presented in binary format (base 2).
  template <typename T, unsigned int const Bits = (8*sizeof(T))>
  struct BinObj: public BitObjHolder<BinObjTag, T, Bits>
    { using BitObjHolder<BinObjTag, T, Bits>::BitObjHolder; };

  struct HexObjTag {}; ///< Tag object for `HexObj`.
  /// Holder for data to be presented in hexadecimal format (base 16).
  template <typename T, unsigned int const Bits = (8*sizeof(T))>
  struct HexObj: public BitObjHolder<HexObjTag, T, Bits>
    { using BitObjHolder<HexObjTag, T, Bits>::BitObjHolder; };

  /**
   * @brief Wrapper to have data printed as hexadecimal dump.
   * @param Atom base type of the dump
   * 
   * This record points to the data to be dumped, and it also include some
   * dumping parameters:
   * * `size`: how many atoms are present in the data
   * * `columns`: how many atoms to print on each line
   * 
   * The data is interpreted as a sequence of `Atom` types.
   * Supported `Atom` types are unsigned integral types (`std::uint8_t`, 
   * `std::uint16_t`, etc.). Each atom is dumped as an `HexObj`.
   */
  template <typename Atom>
  struct HexDumper {
    
    Atom const* const data;
    std::size_t const size;
    unsigned int const columns { 16U };
    
    HexDumper(Atom const* data, std::size_t size, unsigned int columns = 16U)
      : data(data), size(size), columns(columns) {}
    
  }; // HexDumper
  
  
  /**
   * @brief A wrapper padding the dump of its data with zeroes (or `C`).
   * @tparam T type of data to be dumped
   * @tparam Fill (default: `0`) fill character
   * @see `sbnd::ns::util::zeropad()`
   * 
   * The wrapper allows zero-padding of data with a specified field width.
   */
  template <typename T, char Fill = '0'>
  struct ZeroPadder {
    T const data;
    unsigned int const field;
    char const pad;
    
    ZeroPadder(T data, unsigned int field, char pad = Fill)
      : data(data), field(field), pad(pad) {}
    
  }; // struct ZeroPadder
  
  
  /// An object representing `N` characters of value `C`.
  template <unsigned int N, char C = ' '>
  struct Blanks {};
  
  
  /// @}
  // --- END -- Wrapping objects -----------------------------------------------
  
  
  /// Returns a bit mask of type `T` with the four most significant bits set.
  template <typename T>
  constexpr T fourMSBmask();
  
  /// Prints a zero-padded integral `value` into `out`.
  template <typename T>
  void printHex(std::ostream& out, T value);
  
  
  // --- BEGIN -- dump operators -----------------------------------------------
  
  /// Dumps `N` characters of value `C` to `out` stream.
  template <unsigned int N, char C = ' '>
  std::ostream& operator<< (std::ostream& out, Blanks<N, C>);
  
  
  /**
   * @brief Dumps `data` bit by bit into `out` stream.
   * @param out output stream
   * @param data data wrapper with information on how many bits to dump
   * @return the output stream `out`
   * 
   * Parameters of the dump are read from the `data` wrapper.
   * The dump is in format `(Bits) bbb bbbb bbbb ...` (`Bits` is the number
   * of bits, and `b` are bit values, `0` or `1`, the first being the most
   * significant bit).
   */
  template <typename T, unsigned int Bits>
  std::ostream& operator<< (std::ostream& out, BinObj<T, Bits> const& data);
  
  
  /**
   * @brief Dumps `data` nibble by nibble into `out` stream.
   * @param out output stream
   * @param data data wrapper
   * @return the output stream `out`
   * 
   * The value in `data` is printed in hexadecimal format, including all its
   * bits.
   * The STL `std::hex` mode of a stream may be more convenient,
   * but it's sticky:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * using namespace sbnd::ns::util;
   * std::cout << 36 << std::endl; // "36"
   * std::cout << details::HexObj{ 36 } << std::endl; // 00000024
   * std::cout << 36 << std::endl; // "36"
   * std::cout << std::hex << 36 << std::endl; // "24"
   * std::cout << 36 << std::endl; // "24" (std::hex is sticky)
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * Parameters of the dump are read from the `data` wrapper.
   * The dump is in format `(Bits) bbb bbbb bbbb ...` (`Bits` is the number
   * of bits, and `b` are bit values, `0` or `1`, the first being the most
   * significant bit).
   */
  template <typename T>
  std::ostream& operator<< (std::ostream& out, HexObj<T> const& data);
  
  
  /**
   * @brief Dumps data in a hexadecimal table.
   * @param out output stream
   * @param data data wrapper with information on how to dump it
   * @return the output stream `out`
   * 
   * Wrapped data is printed in a table: address of the first `Atom` of data,
   * a separator `|`, a sequence of as many atom values as specified in
   * `data.columns`, in hexadecimal format and zero-padded, and another
   * separator `|`.
   * If there are 6 or more columns, a larger space indentation is inserted
   * between the two central columns.
   * The table is written on a new line, and the line is ended after the table.
   */
  template <typename Atom>
  std::ostream& operator<< (std::ostream& out, HexDumper<Atom> const& data);

  /// Dumps a value padding with `0` via `ZeroPadder` wrapper.
  template <typename T>
  std::ostream& operator<< (std::ostream& out, ZeroPadder<T> const& data);
} // namespace sbnd::ns::util::details


// -----------------------------------------------------------------------------
namespace sbnd::ns::util {
  
  template <typename IOS>
  class FormatFlagsGuard;
  
  
  // --- BEGIN -- Format adapters ----------------------------------------------
  /**
   * @name Format adapters
   * 
   * These functions should be used in stream insertion statements like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << sbnd::ns::util::bin(0xAA) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   */
  /// @{
  
  /**
   * @brief Returns a wrapper to print the specified data in binary format.
   * @tparam T type of datum to be printed
   * @param value the value to be printed
   * @see `sbnd::ns::util::details::operator<< (std::ostream&, sbnd::ns::util::details::BinObj<T, Bits> const&)`   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << sbnd::ns::util::bin(0xAA) << std::endl;
   * std::cout << sbnd::ns::util::bin<char>(0xAA) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print `(32) 0000 0000 0000 0000 0000 0000 1010 1010` on the first line
   * (assuming the representation of `int` is 32 bit) and `(8) 1010 1010` on the
   * second line.
   */
  template <typename T>
  constexpr details::BinObj<T> bin(T value);

  /**
   * @brief Returns a wrapper to print the specified data in binary format.
   * @tparam Bits number of bits to print out of the type `T`
   * @tparam T type of datum to be printed
   * @param value the value to be printed
   * @see `sbnd::ns::util::details::operator<< (std::ostream&, sbnd::ns::util::details::BinObj<T, Bits> const&)`
   * 
   * Only the least significant `Bits` will be printed.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << sbnd::ns::util::bin<10U>(0xAA) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print `(10) 00 1010 1010`.
   */
  template <std::size_t Bits, typename T>
  constexpr details::BinObj<T, Bits> bin(T value);

  /**
   * @brief Returns a wrapper to print the specified data in hex dump format.
   * @tparam Atom type of basic element of the data
   * @param data pointer to the data to be printed
   * @param size number of elements to be printed
   * @param columns (default: `16`) atoms per output line
   * @see `sbnd::ns::util::details::operator<< (std::ostream&, sbnd::ns::util::details::HexDumper<Atom> const&)`
   * 
   * 
   * 
   * Only the least significant `Bits` will be printed.
   * 
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * const char data[] = "012345";
   * std::cout << sbnd::ns::util::hexdump(data, 7U, 8U) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print 7 bytes from `data`, using a 8 column format, with an output
   * similar to:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 0X81234560 | 30 31 32 33  34 35 00    |
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * while
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * const std::uint16_t powers[]
   *   = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096 };
   * std::cout << sbnd::ns::util::hexdump(powers, 13U, 8U) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print 13 values from `data`, using a 8 column format, with an output
   * similar to:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 0X81234570 | 0001 0002 0004 0008  0010 0020 0040 0080 |
   * 0X81234580 | 0100 0200 0400 0800  1000                |
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename Atom>
  details::HexDumper<Atom> hexdump
    (Atom const* data, std::size_t size, unsigned int columns = 16U);

  /**
   * @brief Returns a wrapper to print the specified data with a field width
   * @tparam T type of data to be printed
   * @param data value to be printed
   * @param field number of characters to use
   * @param pad (default: `0`) filling character
   * @return object to be inserted into a stream
   * 
   * The specified `value` is printed right-padded into a space at least `field`
   * character wide, using `pad` as filling character on the left of `value`.
   * 
   * C++ STL I/O is used to produce the output.
   * Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << sbnd::ns::util::zeropad(79, 4) << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print `0079` while
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::cout << std::hex << sbnd::ns::util::zeropad(79, 4, '*') << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will print `**4F`.
   */
  template <typename T>
  details::ZeroPadder<T> zeropad(T data, unsigned int field, char pad = '0');

  /// @}
  // --- END -- Format adapters ------------------------------------------------

} // namespace sbnd::ns::util


// -----------------------------------------------------------------------------
/**
 * @brief Saves some status of the specified stream object, and restores it.
 * @tparam IOS type of stream object, with flag interface like `std::ostream`
 * 
 * This object uses the RIIA pattern to read some status from a stream, and then
 * restore it on the destruction of the object itself.
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::cout << std::hex << 26 << std::endl; // prints "1a"
 * std::cout << 28 << std::endl; // prints "1c" (hexadecimal format sticks)
 * {
 *   sbnd::ns::util::FormatFlagsGuard sg(std::cout); // saves std::cout status
 *   std::cout << std::oct << 26 << std::endl; // prints "32"
 *   std::cout << 28 << std::endl; // prints "34" (octal format sticks)
 *   // destruction of sg happens here, std::cout flags are restored
 * }
 * std::cout << 28 << std::endl; // prints "1c" (hexadecimal format restored)
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * Saved settings currently include:
 * * all format flags (`std::ios_base::flags()`)
 * * fill character (`std::basic_ios::fill()`)
 * 
 */
template <typename IOS>
class sbnd::ns::util::FormatFlagsGuard {
  static_assert(!std::is_const_v<IOS>);
  
  IOS& ios;
  typename IOS::fmtflags const flags;
  typename IOS::char_type const fill;
  
    public:
  FormatFlagsGuard(IOS& ios)
    : ios{ios}, flags{ios.flags()}, fill{ios.fill()} {}
  
  ~FormatFlagsGuard() { restore(); }
  
  void restore() { ios.flags(flags); ios.fill(fill); }
  
}; // struct sbnd::ns::util::FormatFlagsGuard


// -----------------------------------------------------------------------------
template <typename T>
void sbnd::ns::util::details::printHex(std::ostream& out, T value) {
  static_assert(std::is_integral_v<T>, "Only integral types are supported.");
  
  // we actually store the null character at the end; it does not really matter
  constexpr const char HexChars[17U] = "0123456789ABCDEF";
  
  // print nibble by nibble, no spaces, starting from the most significant
  std::size_t nibblesLeft = sizeof(value) * 2U;
  while (nibblesLeft--) {
    std::size_t const digit = (value >> (nibblesLeft * 4U)) & 0xF;
    out << HexChars[digit];
  } // while
} // sbnd::ns::util::details::printHex()


// -----------------------------------------------------------------------------
template <unsigned int N, char C >
std::ostream& sbnd::ns::util::details::operator<<
  (std::ostream& out, Blanks<N, C>)
{
  for (unsigned int i = 0U; i < N; ++i) out.put(C);
  return out;
} // sbnd::ns::util::details::operator<< (sbnd::ns::util::details::Blanks)


// -----------------------------------------------------------------------------
template <typename T, unsigned int Bits>
std::ostream& sbnd::ns::util::details::operator<<
  (std::ostream& out, BinObj<T, Bits> const& data)
{
  static_assert(std::is_integral_v<T>);
  static_assert(Bits > 0U);
  
  unsigned int bitsLeft = data.bits;
  T mask = T{ 1 } << (bitsLeft - 1);
  out << "(" << bitsLeft << ") ";
  while (mask) {
    out << ((data.data & mask)? '1': '0');
    mask >>= 1;
    if (--bitsLeft == 0) break;
    if ((bitsLeft & 0x03) == 0x00) out << ' ';
  } // while
  return out;
} // sbnd::ns::util::details::operator<< (sbnd::ns::util::details::BinObj)


// -----------------------------------------------------------------------------
template <typename T>
std::ostream& sbnd::ns::util::details::operator<<
  (std::ostream& out, HexObj<T> const& data)
  { printHex(out, data.data); return out; }


// -----------------------------------------------------------------------------
template <typename Atom>
std::ostream& sbnd::ns::util::details::operator<<
  (std::ostream& out, HexDumper<Atom> const& data)
{
  
  static constexpr std::size_t AtomChars = sizeof(Atom) * 2;
  static constexpr Blanks<AtomChars> BlankAtom;
  
  auto const printAtoms = [&out]
    (Atom const* ptr, Atom const* const ptrend, std::ptrdiff_t columns)
    {
      /*
      FormatFlagsGuard const outGuard { out };
      out.setf(std::ios_base::hex, std::ios_base::basefield);
      out.unsetf(std::ios_base::showbase);
      out.fill('0');
      */
      Atom const* cend = ptr + columns;
      while (ptr != cend) {
        out << ' ';
        if (ptr < ptrend) {
          // out << std::setw(AtomChars) << *ptr;
          printHex(out, *ptr);
        }
        else              out << BlankAtom;
        ++ptr;
      } // while
      return ptr;
    };
  
  Atom const* ptr = data.data;
  Atom const* const ptrend = ptr + data.size;
  
  FormatFlagsGuard const outGuard { out };
  out.fill('0');
  
  auto const halfColumns = data.columns / 2;
  while (ptr < ptrend) {
    
    out << "\n" << std::setw(8) << ((void*) ptr) << " |";
    
    ptr = printAtoms(ptr, ptrend, data.columns - halfColumns);
    if (data.columns >= 6U) out << ' ';
    ptr = printAtoms(ptr, ptrend, halfColumns);
    out << " |";
    
  } // while
  
  out << '\n';
  
  return out;
} // operator<< (HexDumper)


// -----------------------------------------------------------------------------
template <typename T>
std::ostream& sbnd::ns::util::details::operator<<
  (std::ostream& out, ZeroPadder<T> const& data)
{
  FormatFlagsGuard ffg { out };
  out.fill(data.pad);
  return out << std::setw(data.field) << data.data;
} // sbnd::ns::util::details::operator<< (ZeroPadder)



// -----------------------------------------------------------------------------
template <typename T>
constexpr auto sbnd::ns::util::bin(T value) -> details::BinObj<T>
  { return details::BinObj<T>{ std::move(value) }; }

template <std::size_t Bits, typename T>
constexpr auto sbnd::ns::util::bin(T value) -> details::BinObj<T, Bits>
  { return details::BinObj<T, Bits>{ std::move(value) }; }


// -----------------------------------------------------------------------------
template <typename Atom>
auto sbnd::ns::util::hexdump
  (Atom const* data, std::size_t size, unsigned int columns /* = 16U */)
  -> details::HexDumper<Atom>
  { return details::HexDumper<Atom>{ data, size, columns }; }


// -----------------------------------------------------------------------------
template <typename T>
auto sbnd::ns::util::zeropad(T data, unsigned int field, char pad /* = '0' */)
  -> details::ZeroPadder<T>
  { return details::ZeroPadder<T>{ data, field, pad }; }


// -----------------------------------------------------------------------------

#endif // SBNDCODE_DECODERS_DECODERTOOLS_BINARYDUMPUTILS_H
