/**
 * @file   sbndcode/Decoders/DecoderTools/details/KeyValuesData.h
 * @brief  Simple parsed data format.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * 
 * This library is header only.
 */

#ifndef SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_KEYVALUESDATA_H
#define SBNDCODE_DECODERS_DECODERTOOLS_DETAILS_KEYVALUESDATA_H


// C++ standard libraries
#include <iosfwd> // std::ostream
#include <vector>
#include <string_view>
#include <string>
#include <optional>
#include <stdexcept> // std::runtime_error
#include <utility> // std::move()
#include <charconv> // std::from_chars()
#include <type_traits> // std::is_floating_point_v, std::enable_if_t
#include <cstddef> // std::size_t


// -----------------------------------------------------------------------------
namespace sbnd {
  
  namespace details {
    template <typename T, typename Enable = void>
    struct KeyValuesConverter;
    template <typename T, unsigned int Base = 10U>
    struct BaseWrapper;
  }
  
  struct KeyValuesData;
  
  std::ostream& operator<< (std::ostream& out, KeyValuesData const& data);
  
} // namespace sbnd


// -----------------------------------------------------------------------------
/**
 * @class sbnd::KeyValuesData
 * @brief Collection of items with key/values structure.
 * 
 * This class collects `Item` objects, each being a string key and a sequence
 * of zero or more values. An specific item can be accessed by its key
 * (`findItem()`, `getItem()`) or all items may be iterated through (`items()`).
 * 
 * 
 * Value type conversions
 * -----------------------
 * 
 * The `Item` objects in this class contain unparsed strings. Each `Item` has a
 * key and a sequence of values. The values can be queried as strings or as
 * other data types. Conversions are performed by
 * `sbnd::details::KeyValuesConverter`, which can be specialized with
 * the needed conversion logic. Only one type of conversion to any given type
 * is supported. Alternative conversions may be achieved using type wrappers
 * (e.g. specializing for a `CaseInsensitive<S>` object that contains a string
 * of type `S` and reading/converting into the string the result of the
 * conversion).
 * 
 * Each converter object specialization for a type `T` should support a call
 * with argument `std::string` returning a `std::optional<T>`.
 * 
 * 
 * Initialization and updates
 * ---------------------------
 * 
 * A `KeyValuesData` object always starts empty (implicit default constructor).
 * A new item is also always created empty (`makeItem()`, `makeOrFetchItem()`)
 * and with a set key.
 * 
 * After an empty item (i.e. an item with a key but no values) is created,
 * values can be added using the `Item` subclass interface (`addValue()`).
 * Item values and keys can be modified by changing `Item` data members
 * directly. The item to be modified is immediately returned by `makeItem()`;
 * afterwards, an existing item may be still retrieved for changes with
 * `makeOrFetchItem()` and `findItem()`.
 * 
 * This object will refuse to create a new item with the same key as an existing
 * one. However, the key of the item may be changed to any value after
 * `makeItem()` is called, although the interface does not encourage that.
 * If this introduces a duplicate key, the query functions will systematically
 * retrieve only one of the items with the repeated key (which one and whether
 * always the same one are undefined).
 * 
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * sbnd::KeyValuesData data;
 * data.makeItem("TriggerType").addValue("S5");
 * data.makeItem("Triggers");
 * data.makeItem("TriggerWindows").addValue("0C0B");
 * data.makeItem("TPChits")
 *   .addValue("12").addValue("130").addValue("0").addValue("0");
 * data.makeItem("TPChitTimes")
 *   .addValue("3").addValue("-1.1").addValue("-0.3").addValue("0.1");
 * data.makeItem("PMThits").addValue("8");
 * data.makeOrFetchItem("TPChits").addValue("115");
 * data.makeOrFetchItem("TPChitTimes").addValue("-0.7");
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * creates a `data` object with four items, a `TriggerType` one with one value,
 * a `Triggers` one with no values, a `TriggerWindows` with a single value
 * (expressed as a hexadecimal number), a `TPChits` one with four values,
 * a `TPChitTimes` with four values (the first meant to be the number of
 * remaining ones) and a `PMThits` with one. Finally, it adds one value to
 * `TPChits` and one to `TPChitTimes`, which will both have five afterwards.
 * 
 * 
 * Query
 * ------
 * 
 * The interface is quite terse.
 * 
 * General query methods reports whether there is no item in the object
 * (`empty()`) and how many items there are (`size()`).
 * 
 * A item with a known key can be retrieved (`getItem()`, `findItem()`), or
 * its existence may be tested (`hasItem()`).
 * 
 * Finally, all items can be iterated (`items()`). In this case, the items
 * are presented in the creation order.
 * 
 * The `Item` interface is documented in that subclass.
 * 
 * Example using the `data` object from the previous example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * std::string triggerType = data.getItem("TriggerType").values()[0];
 * std::vector<int> triggers = data.getItem("Triggers").getVector<int>();
 * std::uint32_t triggerWindowBits
 *  = data.getItem("TriggerWindows").getNumber<std::uint32_t>(0, 16); // base 16
 * std::vector<int> TPChits = data.getItem("TPChits").getVector<int>();
 * std::vector<float> TPCtimes
 *  = data.getItem("TPChitTimes").getSizedVector<float>();
 * std::vector<int> CRThits;
 * if (auto const* item = data.findItem("CRThits"))
 *   CRThits = item->getVector<int>();
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Type conversion customization
 * ------------------------------
 *
 * The values in a `Item` object can be queried as strings or as other data
 * types using the `GetAs()` interface. Conversions are performed by
 * `sbnd::details::KeyValuesConverter`, which can be specialized with
 * the needed conversion logic. Only one type of conversion to any given type
 * is supported.
 * 
 * Each converter object specialization for a type `T` should support a call
 * with argument `std::string` returning a `std::optional<T>`.
 * 
 */
struct sbnd::KeyValuesData {
  
  // --- BEGIN --- Exception definition ----------------------------------------
  /// @name Exceptions
  /// @{
  
  struct Error;
  struct ErrorOnKey;
  struct DuplicateKey;
  struct ConversionFailed;
  struct ItemNotFound;
  struct ValueNotAvailable;
  struct MissingSize;
  struct WrongSize;
  
  /// @}
  // --- END ----- Exception definition ----------------------------------------
  
  /**
   * @brief Representation of a single item of data: a key and several values.
   * 
   * Values can be added directly accessing the `values` data member, or
   * with `addValue()` method call.
   * 
   * Access to the values happens by index, with `nValues()` indices starting
   * from `0` available. Direct access to `values()` is expected, and additional
   * shortcuts are available:
   * * `getAs()`, `getOptionalAs()` to convert to an arbitrary type; the
   *     standard conversion uses `from_chars()`, and specialization is
   *     possible
   * * `getNumber()`, `getOptionalNumber()` to convert to a number
   * * `getVector()` to convert to a vector of values (each like in `getAs()`)
   * * `getSizedVector()` to convert to a vector of values; the first value is
   *   the size of the actual values.
   */
  struct Item: public std::pair<std::string, std::vector<std::string>> {
    
    using pair_t = std::pair<std::string, std::vector<std::string>>;
    
    /**
     * @brief Alias for numerical base conversion.
     * @tparam T type of the number to be converted
     * 
     * Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const i16
     *   = item.getAs<int>(0U, sbnd::KeyValuesData::Item::UseBase<int>{ 16 });
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * converts the first value of `item` into an integer (`int`) with base 16.
     */
    template <typename T>
    using UseBase
      = sbnd::details::KeyValuesConverter<sbnd::details::BaseWrapper<T>>;
    
    
    // --- BEGIN -- Constructors -----------------------------------------------
    
    /// Constructs a new item assigning it a key (which should not be changed).
    Item(std::string key): pair_t{ std::move(key), {} } {}
    
    // --- END ---- Constructors -----------------------------------------------
    
    
    // --- BEGIN -- Setting ----------------------------------------------------
    /// @name Setting interface
    /// @{
    
    //@{
    /// Appends a string value to the list of values.
    /// @return this same object (allows queueing calls in the same statement)
    Item& addValue(std::string value)
      { values().push_back(std::move(value)); return *this; }
    Item& addValue(std::string_view value)
      { return addValue(std::string{ value }); }
    //@}
    
    /// Appends a sequence of values to this key.
    template <typename BIter, typename EIter>
    Item& addValues(BIter begin, EIter end)
      { while (begin != end) addValue(*begin++); return *this; }
    
    //@{
    /// Removes all the values.
    void clear() { values().clear(); }
    //@}
    
    /// @}
    // --- END ---- Setting ----------------------------------------------------
    
    
    // --- BEGIN -- Direct access ----------------------------------------------
    /// @name Direct access to key and values
    /// @{
    
    //@{
    /// Returns the key of the item.
    std::string const& key() const noexcept { return first; }
    //@}
    
    //@{
    /// Returns all item values, as strings.
    std::vector<std::string>& values() noexcept { return second; }
    std::vector<std::string> const& values() const noexcept { return second; }
    //@}
    
    //@{
    /// Returns the value at `index` (unchecked) with no conversion.
    std::string const& value(std::size_t index = 0) const noexcept
      { return values()[index]; }
    //@}
    
    //@{
    /// Returns the value at `index` with no conversion, no value if not present.
    std::optional<std::string> optionalValue(std::size_t index) const noexcept;
    //@}
    
    //@{

    /// Returns the number of values currently present.
    std::size_t nValues() const noexcept { return values().size(); }
    
    /// @}
    // --- END ---- Direct access ----------------------------------------------
    
    
    // --- BEGIN -- Query ------------------------------------------------------
    /// @name Query interface
    /// @{
    
    //@{
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param index the index of the requested value
     * @param converter a functor for conversion of the value into type `T`
     * @return the requested value as an object of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     */
    template <typename T, typename Conv>
    T getAs(std::size_t index, Conv converter) const;
    //@}
    
    //@{
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @param index the index of the requested value
     * @return the requested value as an object of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via an helper class
     * `sbnd::details::KeyValuesConverter` which can be specialized if needed,
     * and that uses `from_chars()` for conversion.
     */
    template <typename T>
    T getAs(std::size_t index) const;
    //@}
    
    
    //@{
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam IgnoreFormatErrors (default: `true`) how to treat conversion
     *                            errors
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param index the index of the requested value
     * @param converter a functor for conversion of the value into type `T`
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     * 
     * If no value is available for the specified `index`, an empty optional
     * is returned.
     * 
     * An exception is thrown on conversion failures unless `IgnoreFormatErrors`
     * is `true`, in which case an empty optional is also returned.
     */
    template <typename T, bool IgnoreFormatErrors = false, typename Conv>
    std::optional<T> getOptionalAs(std::size_t index, Conv converter) const;
    //@}
    
    //@{
    /**
     * @brief Returns the requested value, converted into type `T`
     * @tparam T type to convert the value into
     * @tparam IgnoreFormatErrors (default: `true`) how to treat conversion
     *                            errors
     * @param index the index of the requested value
     * @param ignoreFormatErrors (default: `false`) ignore conversion errors
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * Conversion is performed via `converter` object, functor taking a string
     * and returning an object of type `std::optional<T>`. The functor can
     * decline the conversion by returning an empty `std::optional`, or directly
     * throw an exception on error.
     * 
     * If no value is available for the specified `index`, an empty optional
     * is returned.
     * 
     * An exception is thrown on conversion failures unless `IgnoreFormatErrors`
     * is `true`, in which case an empty optional is also returned.
     */
    template <typename T, bool IgnoreFormatErrors = false>
    std::optional<T> getOptionalAs(std::size_t index) const;
    //@}
    
    
    //@{
    /**
     * @brief Returns the requested value, converted into a number of type `T`
     * @tparam T type of number to convert the value into
     * @param index the index of the requested value
     * @param base (default: `10`) numerical base of the input number
     * @return the requested value as a number of type `T`
     * @throw ValueNotAvailable if no value is available with that index
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * See `getAs()` for details.
     * 
     * Note that the number must have no base prefix (e.g. `"F5"` for
     * hexadecimal rather than `"0xF5"`).
     */
    template <typename T>
    T getNumber(std::size_t index, unsigned int base) const;
    template <typename T>
    T getNumber(std::size_t index) const;
    //@}
    
    // TODO do we need to propagate the IgnoreFormatErrors option?
    //@{
    /**
     * @brief Returns the requested value, converted into a number of type `T`
     * @tparam T type of number to convert the value into
     * @param index the index of the requested value
     * @param base (default: `10`) numerical base of the input number
     * @return the requested value, or an empty optional on failure
     * @throw ConversionFailed if the value could not be converted to type `T`
     * 
     * See `getOptionalAs()` for details.
     * 
     * Note that the number must have no base prefix (e.g. `"F5"` for
     * hexadecimal rather than `"0xF5"`).
     */
    template <typename T>
    std::optional<T> getOptionalNumber
      (std::size_t index, unsigned int base) const;
    template <typename T>
    std::optional<T> getOptionalNumber(std::size_t index) const;
    //@}
    
    
    //@{
    /**
     * @brief Returns all the values, each converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param converter a functor for conversion of the value into type `T`
     * @return a vector with all the converted values
     * @throw ConversionFailed if any value could not be converted to type `T`
     * 
     * Conversion of each element is performed by `getAs<T, Conv>()`.
     * 
     * An exception is thrown on any conversion failure.
     */
    template <typename T, typename Conv = details::KeyValuesConverter<T>>
    std::vector<T> getVector(Conv converter = {}) const;
    //@}
    
    //@{
    /**
     * @brief Returns all the values, each converted into type `T`
     * @tparam T type to convert the value into
     * @tparam Conv type of a functor for conversion of the value into type `T`
     * @param converter a functor for conversion of the value into type `T`
     * @return a vector with all the converted values
     * @throw MissingSize on any error converting the first value to a size
     * @throw WrongSize if the actual number of values does not match the size
     * @throw ConversionFailed if any value could not be converted to type `T`
     * 
     * The first value (mandatory) is converted to represent the size of the
     * vector. That is used as verification when converting all the other
     * elements: if there is the wrong number of elements, an exception is
     * thrown.
     * 
     * Conversion of each element is performed by `getAs<T, Conv>()`.
     * 
     * An exception is also thrown on conversion failure of any of the values.
     */
    template <typename T, typename Conv = details::KeyValuesConverter<T>>
    std::vector<T> getSizedVector(Conv converter = Conv{}) const;
    //@}
    
    /// @}
    // --- END ---- Query ------------------------------------------------------
    
    
    /// Lexicographic order by key (case-sensitive).
    bool operator< (Item const& other) const noexcept
      { return key() < other.key(); }
    
    
      private:
    
    // implementation detail
    template <typename T, typename Iter, typename Conv>
    std::vector<T> convertVector(Iter begin, Iter end, Conv converter) const;

    /// Conversion functions.
    template <typename T>
    std::optional<T> convertStringInto(std::string const& valueStr) const
      { return sbnd::details::KeyValuesConverter<T>{}(valueStr); }
    
  }; // struct Item
  
  
  // --- BEGIN -- Setter interface ---------------------------------------------
  /// @name Setter interface
  /// @{
  
  /// @brief Creates and registers a new item with the specified `key`.
  /// @return the newly created item for modifications
  Item& makeItem(std::string key);
  
  /// @brief Creates or retrieves an item with the specified `key`.
  /// @return the newly created or existing item for modifications
  Item& makeOrFetchItem(std::string const& key);
  
  /// Returns the item with specified `key`, `nullptr` if none.
  Item* findItem(std::string const& key) noexcept;
  
  /// @}
  // --- END ---- Setter interface ---------------------------------------------
  
  
  // --- BEGIN -- Query interface ----------------------------------------------
  /// @name Query interface
  /// @{
  
  /// Returns the item with specified `key`, `nullptr` if none.
  Item const* findItem(std::string const& key) const noexcept;
  
  /// Returns the item with specified `key`, throws `std::out_of_range` if none.
  Item const& getItem(std::string const& key) const;
  
  /// Returns whether an item with the specified key is present.
  bool hasItem(std::string const& key) const noexcept;
  
  /// Returns whether there is no item in data.
  bool empty() const noexcept;
  
  /// Returns the number of items in the data.
  std::size_t size() const noexcept;
  
  /// Returns a forward-iterable list of references to items.
  decltype(auto) items() const noexcept;
  
  /// @}
  // --- END ---- Query interface ----------------------------------------------
  
    private:
  
  std::vector<Item> fItems; ///< Collection of data items.
  
  
  /// Creates, registers and return a new item (assumed not to exist yet).
  Item& makeItemImpl(std::string key);
  
}; // sbnd::KeyValuesData


// -----------------------------------------------------------------------------
namespace sbnd {
  std::ostream& operator<<
    (std::ostream& out, KeyValuesData::Item const& data);
} // namespace sbnd


// -----------------------------------------------------------------------------
// ---  Exception class definitions
// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::Error: public std::runtime_error {
  
  Error(std::string msg): std::runtime_error(std::move(msg)) {}
  
}; // sbnd::KeyValuesData::Error()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::ErrorOnKey: public Error {
  
  ErrorOnKey(std::string const& key, std::string const& msg)
    : Error("Key '" + key + "': " + msg) {}
  
}; // sbnd::KeyValuesData::ErrorOnKey()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::DuplicateKey: public Error {
  
  DuplicateKey(std::string const& msg)
    : Error("KeyValuesData::DuplicateKey: '" + msg + '\'')
    {}
  
}; // sbnd::KeyValuesData::DuplicateKey()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::ConversionFailed: public ErrorOnKey {
  
  ConversionFailed(
    std::string const& key, std::string const& s, std::string const& tname = ""
    )
    : ErrorOnKey{ key,
      "conversion of '" + s + "'"
       + (tname.empty()? "": (" to type '" + tname + "'")) + " failed"
      }
    {}
  
  template <typename T>
  static ConversionFailed makeFor(std::string const& key, std::string const& s)
    { return { key, s, typeid(T).name() }; }
  
  template <typename T>
  static ConversionFailed makeFor
    (std::string const& key, std::size_t index, std::string const& s)
    { return makeFor<T>(key + "[" + std::to_string(index) + "]", s); }
  
}; // sbnd::KeyValuesData::ConversionFailed()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::ItemNotFound: public ErrorOnKey {
  
  ItemNotFound(std::string const& key): ErrorOnKey(key, "key not found") {}
  
}; // sbnd::KeyValuesData::ItemNotFound()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::ValueNotAvailable: public ErrorOnKey {
  
  ValueNotAvailable(std::string const& key, std::size_t index)
    : ErrorOnKey(key, "item value #" + std::to_string(index) + " not available")
    {}
  
}; // sbnd::KeyValuesData::ValueNotAvailable()


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::MissingSize: public ErrorOnKey {
  
  MissingSize(std::string const& key)
    : ErrorOnKey
      (key, "is required to have a size as first value, but it has no values")
    {}
  
  MissingSize(std::string const& key, std::string const& valueStr)
    : ErrorOnKey(
      key,
      " first value '" + valueStr + "' can't be converted into a vector size"
      )
    {}
  
}; // sbnd::KeyValuesData::MissingSize


// -----------------------------------------------------------------------------
struct sbnd::KeyValuesData::WrongSize: public ErrorOnKey {
  
  WrongSize(std::string const& key, std::size_t expected, std::size_t actual)
    : ErrorOnKey(key,
      std::to_string(expected) + " values (except the size) were expected, "
      + std::to_string(actual) + " are present instead"
      )
    {}
  
}; // sbnd::KeyValuesData::WrongSize


// -----------------------------------------------------------------------------
// --- Template implementation
// -----------------------------------------------------------------------------
namespace sbnd::details {
  
  template <typename T, unsigned int Base /* = 10U */>
  struct BaseWrapper {
    T value;
    constexpr operator T() const noexcept { return value; }
  }; // BaseWrapper
  
} // namespace sbnd::details


// -----------------------------------------------------------------------------
// ---  sbnd::KeyValuesData::Item
// -----------------------------------------------------------------------------
inline std::optional<std::string> sbnd::KeyValuesData::Item::optionalValue
  (std::size_t index) const noexcept
{
  return (index < nValues())? std::optional{ values()[index] }: std::nullopt;
}


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
T sbnd::KeyValuesData::Item::getAs
  (std::size_t index, Conv converter) const
{
  
  if (index >= values().size()) throw ValueNotAvailable(key(), index);
  
  auto const& valueStr = values()[index];
  auto const number = converter(valueStr);
  return number? *number: throw ConversionFailed::makeFor<T>(key(), valueStr);
  
} // sbnd::KeyValuesData::Item::getAs<>()


// -----------------------------------------------------------------------------
template <typename T>
T sbnd::KeyValuesData::Item::getAs(std::size_t index) const
  { return getAs<T>(index, details::KeyValuesConverter<T>{}); }


// -----------------------------------------------------------------------------
template <typename T, bool IgnoreFormatErrors /* = false */, typename Conv>
std::optional<T> sbnd::KeyValuesData::Item::getOptionalAs
  (std::size_t index, Conv converter) const
{
  if (index < values().size()) return std::nullopt;

  auto const& valueStr = values()[index];
  auto const number = converter(valueStr);
  return (number || IgnoreFormatErrors)
    ? number: throw ConversionFailed::makeFor<T>(key(), valueStr);

} // sbnd::KeyValuesData::Item::getOptionalAs()


// -----------------------------------------------------------------------------
template <typename T, bool IgnoreFormatErrors /* = false */>
std::optional<T> sbnd::KeyValuesData::Item::getOptionalAs
  (std::size_t index) const
{
  return getOptionalAs<T, IgnoreFormatErrors>
    (index, details::KeyValuesConverter<T>{});
}


// -----------------------------------------------------------------------------
template <typename T>
T sbnd::KeyValuesData::Item::getNumber
  (std::size_t index, unsigned int base) const
  { return getAs<T>(index, UseBase<T>{ base }); }


// -----------------------------------------------------------------------------
template <typename T>
T sbnd::KeyValuesData::Item::getNumber(std::size_t index) const
  { return getAs<T>(index, details::KeyValuesConverter<T>{}); }


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> sbnd::KeyValuesData::Item::getOptionalNumber
  (std::size_t index, unsigned int base) const
  { return getOptionalAs<T>(index, UseBase<T>{ base }); }


// -----------------------------------------------------------------------------
template <typename T>
std::optional<T> sbnd::KeyValuesData::Item::getOptionalNumber
  (std::size_t index) const
  { return getOptionalAs<T>(index, details::KeyValuesConverter<T>{}); }


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
std::vector<T> sbnd::KeyValuesData::Item::getVector
  (Conv converter /* = {} */) const
{
  return
    convertVector<T>(values().begin(), values().end(), std::move(converter));
}


// -----------------------------------------------------------------------------
template <typename T, typename Conv>
std::vector<T> sbnd::KeyValuesData::Item::getSizedVector
  (Conv converter /* = {} */) const
{
  
  if (values().empty()) throw MissingSize(key());
  
  std::size_t const n = getNumber<std::size_t>(0U);
  if (n != values().size() - 1)
    throw WrongSize(key(), n, values().size() - 1);
  
  return convertVector<T>
    (std::next(values().begin()), values().end(), std::move(converter));
} // sbnd::KeyValuesData::Item::getSizedVector()


// -----------------------------------------------------------------------------
template <typename T, typename Iter, typename Conv>
std::vector<T> sbnd::KeyValuesData::Item::convertVector
  (Iter begin, Iter end, Conv converter) const
{
  std::vector<T> data;
  data.reserve(std::distance(begin, end));
  Iter it = begin;
  while (it != end) {
    std::string const& valueStr = *it;
    if (std::optional const number = converter(valueStr))
      data.push_back(*number);
    else {
      throw ConversionFailed::makeFor<T>
        (key(), std::distance(begin, it), valueStr);
    }
    ++it;
  } // while
  return data;
} // sbnd::KeyValuesData::Item::convertVector()


// -----------------------------------------------------------------------------
inline decltype(auto) sbnd::KeyValuesData::items() const noexcept
  { return fItems; }


// -----------------------------------------------------------------------------
// ---  sbnd::details::KeyValuesConverter
// -----------------------------------------------------------------------------
template <typename T, typename /* = void */>
struct sbnd::details::KeyValuesConverter {
  /*
   * The "generic" converter supports arithmetic values and stuff that can be
   * converted to string. For any other type and conversions, specialization
   * is needed.
   * 
   * NOTE: until we adopt C++17-compliant compilers (or better, libraries),
   *       also floating point numbers require specialization.
   */
  
  
  //@{
  /// Convert a string `s` into a type `T`;
  /// may return `std::nullopt` on "non-fatal" failure.
  std::optional<T> operator() (std::string const& s) const
    { return convert(s); }
  
  static std::optional<T> convert(std::string const& s)
    {
      if constexpr (std::is_arithmetic_v<T>) {
        T number {}; // useless initialization to avoid GCC complains
        char const *b = s.data(), *e = b + s.length();
        return (std::from_chars(b, e, number).ptr == e)
          ? std::make_optional(number): std::nullopt;
      }
      else if constexpr(std::is_constructible_v<T, std::string>){
        return std::make_optional(T{ s });
      }
      else return std::nullopt;
    } // convert()
  
  //@}
}; // sbnd::details::KeyValuesConverter<>


/// Specialization for conversions with a numeric base.
template <typename T, unsigned int Base, typename Enable>
struct sbnd::details::KeyValuesConverter
  <sbnd::details::BaseWrapper<T, Base>, Enable>
{
  // NOTE: this specialization returns `T`, not BaseWrapper<T>.
  
  unsigned int base { Base }; ///< The numerical base this converter uses.
  
  //@{
  /// Convert a string `s` into a numerical type `T` from the specified `base`;
  /// may return `std::nullopt` on "non-fatal" failure.
  std::optional<T> operator() (std::string const& s) const
    { return convert(s); }
  
  std::optional<T> convert(std::string const& s) const
    {
      T number {}; // useless initialization to avoid GCC complains
      char const *b = s.data(), *e = b + s.length();
      return (std::from_chars(b, e, number, base).ptr == e)
        ? std::make_optional(number): std::nullopt;
    } // convert()
  
  //@}
}; // sbnd::details::KeyValuesConverter<>


// ---  BEGIN  --- WORKAROUND --------------------------------------------------
/*
 * Compilers like GCC 9.3.0 and Clang 7.0 are not fully C++17 compliant yet.
 * `std::from_chars()` is not provided for floating points types.
 * The first compiler versions supporting them are GCC 11.1 and Clang 12.0.1.
 * I am providing a specialization to cover from that, until the compilers
 * are updated.
 */
#ifdef __GNUC__
#  if (__GNUC__ >= 11)
     // GCC 11 should support std::from_chars() for floating point types;
     // remove this #if branch, and if Clang's is also already removed,
     // remove the workaround too
#    error "Redundant workaround on std::from_chars() for GCC"
#  else
#    define SBNDCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT
#  endif // __GNUC__
#endif // __GNUC__

#ifdef __clang_major__
#  if (__clang_major__ >= 12) || ((__clang_major__ == 12) && ((__clang_minor__ >= 1) || (__clang_patchlevel__ >= 1)))
     // Clang 12.0.1 should support std::from_chars() for floating point types;
     // remove this #if branch, and if GCC's is also already removed,
     // remove the workaround too
#    error "Redundant workaround on std::from_chars() for GCC"
#  else
#    define SBNDCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT
#  endif // __clang_major__
#endif // __clang_major__

#ifdef SBNDCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_NEEDS_FROMCHARS_FLOAT

#include <sstream>

template <typename T>
struct sbnd::details::KeyValuesConverter
  <T, std::enable_if_t<std::is_floating_point_v<T>>>
{
  
  std::optional<T> operator() (std::string const& s) const
    { return convert(s); }
  
  static std::optional<T> convert(std::string const& s)
    {
      T number {}; // useless initialization to avoid GCC complains
      std::istringstream sstr{ s };
      sstr >> number;
      // check that no non-space character is left in the stream
      return (sstr && (sstr >> std::ws).eof())
        ? std::make_optional(number): std::nullopt;
    } // convert()
  
}; // sbnd::details::KeyValuesConverter<floating point>

#endif
// ---  END ------ WORKAROUND --------------------------------------------------

// -----------------------------------------------------------------------------

#endif // SBNDCODE_DECODE_DECODERTOOLS_DETAILS_KEYVALUESDATA_H
