/**
 * @file   sbndcode/Decoders/DecoderTools/details/KeyedCSVparser.h
 * @brief  Simple parser for comma-separated text.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * @see    sbndcode/Decoders/DecoderTools/details/KeyedCSVparser.cxx
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H

// SBND libraries
#include "sbndcode/Decoders/DecoderTools/details/KeyValuesData.h"

// C++ standard libraries
#include <iosfwd> // std::ostream
#include <string_view>
#include <vector>
#include <string>
#include <optional>
#include <regex>
#include <initializer_list>
#include <stdexcept> // std::runtime_error
#include <utility> // std::move(), std::pair
#include <limits>
#include <type_traits> // std::is_constructible_v, std::is_arithmetic_v, ...
#include <charconv> // std::from_chars()
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace sbnd::details { class KeyedCSVparser; }
/**
 * @class sbnd::details::KeyedCSVparser
 * @brief Parser to fill a `KeyValuesData` structure out of a character buffer.
 * 
 * It currently supports only single-line buffer.
 * 
 * The parser operates one "line" at a time, returning a `KeyValuesData` with
 * the values assigned to each detected key. No data type is implied: all
 * elements are treated as strings, either a key or a value.
 * The parser separates the elements according to a separator, strips them of
 * trailing and heading spaces, then it decides whether each element is a value
 * to be assigned to the last key found, or a new key.
 * Keys are elements that have letters in them, values are anything else.
 * This simple (and arguable) criterion can be broken with specific parser
 * configuration: a pattern can be specified that when matched to an element
 * will make it a key; the pattern can also set the number of values that key
 * will require.
 * 
 * For example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * sbnd::details::KeyedCSVparser parser;
 * parser.addPatterns({
 *       { "TriggerType", 1U } // expect one value (even if contains letters)
 *     , { "TriggerWindows", 1U } // expect one value (even if contains letters)
 *     , { "TPChitTimes", sbnd::details::KeyedCSVparser::FixedSize }
 *          // the first value is an integer, count of how many other values
 *   });
 * 
 * sbnd::KeyValuesData data = parser(
 *   "TriggerType, S5, Triggers, TriggerWindows, 0C0B,"
 *   " TPChits, 12, 130, 0, 0, TPChitTimes, 3, -1.1, -0.3, 0.1, PMThits, 8"
 *   );
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * will return `data` with 6 items.
 */
class sbnd::details::KeyedCSVparser {
  
    public:
  
  using ParsedData_t = sbnd::KeyValuesData;
  
  /// Base of all errors by KeyedCSVparser.
  using Error = sbnd::KeyValuesData::Error;
  using ErrorOnKey = sbnd::KeyValuesData::ErrorOnKey;
  struct ParserError; ///< Generic error: base of all errors by KeyedCSVparser.
  struct InvalidFormat; ///< Parsing format is not understood.
  /// Expected number of values is missing.
  using MissingSize = KeyValuesData::MissingSize;
  struct MissingValues; ///< Expected values are missing.
  
  /// Mnemonic size value used in `addPattern()` calls.
  static constexpr unsigned int FixedSize
    = std::numeric_limits<unsigned int>::max();
  /// Mnemonic size value used in `addPattern()` calls.
  static constexpr unsigned int DynamicSize = FixedSize - 1U;
  
  
  /// Constructor: specifies the separator character.
  KeyedCSVparser(char sep = ','): fSep(sep) {}
  
  //@{
  /// Parses the buffer `s` and returns a data structure with the content.
  ParsedData_t parse(std::string_view const& s) const;
  ParsedData_t parse(std::string const& s) const;
  template <typename BIter, typename EIter>
  ParsedData_t parse(BIter b, EIter e) const;
  
  ParsedData_t operator() (std::string_view const& s) const { return parse(s); }
  ParsedData_t operator() (std::string const& s) const { return parse(s); }
  template <typename BIter, typename EIter>
  ParsedData_t operator() (BIter b, EIter e) const { return parse(b, e); }
  //@}
  
  //@{
  /// Parses the buffer `s` and fills `data` with it.
  void parse(std::string_view const& s, ParsedData_t& data) const;
  //@}
  
  /**
   * @name Know patterns
   * 
   * The parser normally treats as a value everything that does not start with a
   * letter.
   * Known patterns may override this behaviour: if a token matches a known
   * pattern, it is considered a key and it is possible to specify the expected
   * number of values.
   * 
   * The number of values can be:
   * * a number: exactly that number of values are required; an exception will
   *   be thrown if not enough tokens are available;
   * * `FixedSize`: the next token must be a non-negative integer specifying how
   *   many other values to add (read this with `Item::getSizedVector()`);
   *   an exception will be thrown if not enough tokens are available;
   * * `DynamicSize`: the standard algorithm is used and values are added as
   *   long as they don't look like keys; the token matching the pattern is
   *   interpreted as a key though.
   * 
   * Patterns are considered in the order they were added.
   */
  /// @{
  
  //@{
  /**
   * @brief Adds a single known pattern.
   * @param pattern the regular expression matching the key for this pattern
   * @param values the number of values for this pattern
   * @return this parser (`addPattern()` calls may be chained)
   */
  KeyedCSVparser& addPattern(std::regex pattern, unsigned int values)
    { fPatterns.emplace_back(std::move(pattern), values); return *this; }
  KeyedCSVparser& addPattern(std::string const& pattern, unsigned int values)
    { return addPattern(std::regex{ pattern }, values); }
  //@}
  
  //@{
  /**
   * @brief Adds known patterns.
   * @param patterns sequence of patterns to be added
   * @return this parser (`addPatterns()` calls may be chained)
   * 
   * Each pattern is a pair key regex/number of values, like in `addPattern()`.
   */
  KeyedCSVparser& addPatterns
    (std::initializer_list<std::pair<std::regex, unsigned int>> patterns);
  KeyedCSVparser& addPatterns
    (std::initializer_list<std::pair<std::string, unsigned int>> patterns);
  //@}
  
  /// @}
  
    private:
  using Buffer_t = std::string_view;
  using SubBuffer_t = std::string_view;
  
  char const fSep = ','; ///< Character used as token separator.
  
  /// List of known patterns for matching keys, and how many values they hold.
  std::vector<std::pair<std::regex, unsigned int>> fPatterns;
  
  /// Returns the length of the next toke, up to the next separator (excluded).
  std::size_t findTokenLength(Buffer_t const& buffer) const noexcept;

  /// Returns the value of the next token, stripped.
  SubBuffer_t peekToken(Buffer_t const& buffer) const noexcept;
  
  /// Extracts the next token from the `buffer` and returns its value, stripped.
  SubBuffer_t extractToken(Buffer_t& buffer) const noexcept;
  
  /// Is content of `buffer` a key (as opposed to a value)?
  bool isKey(SubBuffer_t const& buffer) const noexcept;
  
  
  template <typename String>
  static Buffer_t makeBuffer(String const& s) noexcept;

  static Buffer_t& moveBufferHead(Buffer_t& buffer, std::size_t size) noexcept;
  
  static SubBuffer_t strip(SubBuffer_t s) noexcept;
  static SubBuffer_t stripLeft(SubBuffer_t s) noexcept;
  static SubBuffer_t stripRight(SubBuffer_t s) noexcept;
  static SubBuffer_t stripRightChar(SubBuffer_t s, char c) noexcept;
  
  template <char... Chars>
  static SubBuffer_t stripRightChars(SubBuffer_t s) noexcept;
  
}; // sbnd::details::KeyedCSVparser



// -----------------------------------------------------------------------------
// ---  Exception class definitions
// -----------------------------------------------------------------------------
struct sbnd::details::KeyedCSVparser::ParserError: public Error {
  
  ParserError(std::string msg): Error(std::move(msg)) {}
  
}; // sbnd::details::KeyedCSVparser::ParseError


// -----------------------------------------------------------------------------
struct sbnd::details::KeyedCSVparser::InvalidFormat: public ParserError {
  
  InvalidFormat(std::string const& msg): ParserError("Format error: " + msg) {}
  
}; // sbnd::details::KeyedCSVparser::InvalidFormat


// -----------------------------------------------------------------------------
struct sbnd::details::KeyedCSVparser::MissingValues: public ErrorOnKey {
  
  MissingValues(std::string const& key, unsigned int values)
    : ErrorOnKey(key,
      "data ended while expecting " + std::to_string(values) + " more values"
      )
    {}
  
}; // sbnd::details::KeyedCSVparser::MissingValues


// -----------------------------------------------------------------------------
// --- Inline implementation
// -----------------------------------------------------------------------------
inline auto sbnd::details::KeyedCSVparser::parse
  (std::string_view const& s) const -> ParsedData_t
  { ParsedData_t data; parse(s, data); return data; }


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
template <typename BIter, typename EIter>
auto sbnd::details::KeyedCSVparser::parse(BIter b, EIter e) const
  -> ParsedData_t
  { return parse(std::string_view{ &*b, std::distance(b, e) }); }


// -----------------------------------------------------------------------------


#endif // SBNDCODE_DECODE_DECODERTOOLS_DETAILS_KEYEDCVSPARSER_H
