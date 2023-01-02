/**
 * @file   sbndcode/Decoders/DecoderTools/details/KeyValuesData.cxx
 * @brief  Simple parsed data format.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * @see    sbndcode/Decoders/DecoderTools/details/KeyValuesData.h
 */

// SBND libraries
#include "sbndcode/Decoders/DecoderTools/details/KeyValuesData.h"

// C++ standard libraries
#include <ostream>
#include <utility> // std::move(), std::as_const()


// -----------------------------------------------------------------------------
// ---  sbnd::KeyValuesData
// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::makeItem(std::string key) -> Item& {
  if (hasItem(key)) throw DuplicateKey{ std::move(key) };
  return makeItemImpl(std::move(key));
} // sbnd::KeyValuesData<>::makeItem()


// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::makeOrFetchItem(std::string const& key) -> Item& {
  
  Item* item = findItem(key);
  return item? *item: makeItemImpl(key);
  
} // sbnd::KeyValuesData<>::makeOrFetchItem()


// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::findItem
  (std::string const& key) const noexcept -> Item const*
{
  for (auto const& item: fItems) if (key == item.key()) return &item;
  return nullptr;
} // sbnd::KeyValuesData<>::findItem() const


// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::findItem
  (std::string const& key) noexcept -> Item*
{
  // no violations here: this is a non-const method, with the right to modify
  // object data; and this avoids code duplication.
  return const_cast<Item*>(std::as_const(*this).findItem(key));
} // sbnd::KeyValuesData<>::findItem()


// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::getItem(std::string const& key) const
  -> Item const&
{
  if (auto item = findItem(key); item) return *item;
  throw ItemNotFound(key);
} // sbnd::KeyValuesData<>::getItem()


// -----------------------------------------------------------------------------
bool sbnd::KeyValuesData::hasItem
  (std::string const& key) const noexcept
{
  return findItem(key);
} // sbnd::KeyValuesData<>::hasItem()


// -----------------------------------------------------------------------------
bool sbnd::KeyValuesData::empty() const noexcept
  { return fItems.empty(); }


// -----------------------------------------------------------------------------
std::size_t sbnd::KeyValuesData::size() const noexcept
  { return fItems.size(); }


// -----------------------------------------------------------------------------
auto sbnd::KeyValuesData::makeItemImpl(std::string key) -> Item& {
  fItems.emplace_back(std::move(key));
  return fItems.back();
} // sbnd::KeyValuesData<>::makeItemImpl()


// -----------------------------------------------------------------------------
std::ostream& sbnd::operator<<
  (std::ostream& out, KeyValuesData::Item const& item)
{
  out << "'" << item.key() << "' (" << item.nValues() << ")";
  if (!item.values().empty()) {
    out << ':';
    for (auto const& value: item.values()) out << " '" << value << '\'';
  }
  return out;
} // sbnd::operator<< (KeyValuesData::Item)


// -----------------------------------------------------------------------------
std::ostream& sbnd::operator<< (std::ostream& out, KeyValuesData const& data)
{
  out << data.size() << " items";
  if (!data.empty()) {
    out << ':';
    for (auto const& item: data.items()) out << "\n  " << item;
  } // if has data
  return out << "\n";
} // sbnd::operator<< (KeyValuesData)


// -----------------------------------------------------------------------------
