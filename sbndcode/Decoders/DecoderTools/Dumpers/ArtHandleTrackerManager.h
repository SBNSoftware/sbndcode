/**
 * @file   sbndcode/Decoders/DecoderTools/Dumpers/ArtHandleTrackerManager.h
 * @brief  Tracks handles for cache deletion.
 * @author Afroditi Papadopoulou (apapadopoulou@anl.gov)
 * @date   Jan 28, 2023
 * 
 * This is a header-only library.
 */

#ifndef SBNDCORE_DECODERS_DECODERTOOLS_DUMPERS_ARTHANDLETRACKERMANAGER_H
#define SBNDCORE_DECODERS_DECODERTOOLS_DUMPERS_ARTHANDLETRACKERMANAGER_H


// framework libraries
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "canvas/Utilities/InputTag.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// C/C++ standard libraries
#include <algorithm> // std::count_if()
#include <vector>
#include <memory> // std::unique_ptr<>
#include <any>
#include <utility> // std::forward()
#include <typeinfo>


// -----------------------------------------------------------------------------
namespace util {
  template <typename Event> struct ArtHandleTrackerInterface;
  template <typename Event> class ArtHandleTrackerManager;
  template <typename Event> class LocalArtHandleTrackerManager;
}

// -----------------------------------------------------------------------------
/**
 * @brief Manages handle trackers for an easy call of `removeCachedProduct()`.
 * @tparam Event the type of data product source (the "principal") to serve
 * 
 * This handle manager is designed to simplify the usage of
 * `art::Event::removeCachedProduct()` on multiple data products.
 * 
 * The envisioned usage model in a single-thread module is the following:
 * 1. The manager is a data member of the module (in alternative, it should be
 *    passed down to the places where data products are read, and at a minimum
 *    it needs to register all the handles it is supposed to "manage").
 * 2. The manager can ask the event to read a data product anew, and get an
 *    handle for it. It will register the handle, and return it.
 *     * The manager can also register an existing handle.
 * 3. The manager can deliver a stored handle (but that's a slow process in the
 *    current implementation and it requires the input tag to resolve
 *    ambiguities).
 * 4. Upon request, all handles are asked to remove their cached products.
 *    This request presumably happens at the end of the event-processing
 *    function (`produce()`, `filter()`, `analyze()`).
 *    Note that after an handle has its cached data removed, it's `clear()`'ed.
 * 
 * An handle manager can serve only one _art_ event at a time.
 * All handles come from that event (or are assumed to).
 * 
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * class MyModule: public art::EDAnalyzer {
 *   
 *   util::ArtHandleTrackerManager<art::Event> fDataCacheRemover;
 *   
 *   art::InputTag const fTag;
 *   
 *   // ...
 *   
 *   void analyze(art::Event const & event) override
 *     {
 *       fDataCacheRemover.useEvent(event);
 *       
 *       auto const& handle = fDataCacheRemover.getHandle<MyDataProduct>(fTag);
 *       
 *       auto results = processData(*handle); // ... or whatever
 *       
 *       fDataCacheRemover.removeCachedProducts(); // free caches after use
 *       
 *       results.Write(); // ... or another portion of whatever
 *     }
 *   
 * };
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * Or one can create the handle as preferred, and then register it:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * class MyModule: public art::EDAnalyzer {
 *   
 *   util::ArtHandleTrackerManager<art::Event> fDataCacheRemover;
 *   
 *   art::InputTag const fTag;
 *   
 *   // ...
 *   
 *   void analyze(art::Event const & event) override
 *     {
 *       fDataCacheRemover.useEvent(event);
 *       
 *       auto const& handle = event.getValidHandle<MyDataProduct>(fTag);
 *       
 *       fDataCacheRemover.registerHandle(handle);
 *       
 *       auto results = processData(*handle); // ... or whatever
 *       
 *       fDataCacheRemover.removeCachedProducts(); // free caches after use
 *       
 *       results.Write(); // ... or another portion of whatever
 *     }
 *   
 * };
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * 
 * 
 * Technical notes
 * ================
 * 
 * Support for different types of handles
 * ---------------------------------------
 * 
 * Currently, both `art::Handle` and `art::ValidHandle` (and nothing else)
 * can be created by `ArtHandleTrackerManager`.
 * The code can be extended to support any object that can be
 * passed to `art::Event::removeCachedProduct()` (e.g. `art::ProductID`,
 * if it will ever happen).
 * 
 * The current implementation is expected to be able to recognize handles
 * pointing to the same data product even if they are of different type
 * (e.g. `art::Handle<T>` and `art::ValidHandle<T>`).
 * The manager stores a copy of the first type of handle registered for a given
 * data product; then, if another handle of any type to the same data product
 * is requested or registered, `registerHandle()` will return the same handle
 * in argument, and the `getHandle()` family of functions will get and return a
 * new handle. In both cases, no new registration will happen.
 * 
 * 
 * Multithreading
 * ---------------
 * 
 * The current implementation is inherently not safe for _art_ multithreading.
 * While the inner data structure don't need global state, and the interface
 * can easily be extended to allow for no global state as well, the operations
 * are still performed on all the registered handles at once (meaning, when one
 * thread asks for `removeCachedProduct()` or `forgetAllHandles()`, data from
 * all events is affected).
 * This can be overcome by changing the internal storage (to be reentrant and
 * possibly by event) and a bit of the interface.
 * 
 * The object as it is now can be implemented on top of such object, preserving
 * the current event information and delivering it to the new manager under the
 * hood, with minimal overhead.
 * 
 * If such feature is needed, ask the author (and be prepared to test it).
 * 
 */
template <typename Event>
class util::ArtHandleTrackerManager {
  
    public:
  
  using Event_t = Event; ///< Type of data viewer object to operate with.
  
  /// Configuration record.
  struct Config_t {
    
    /// Name of the output category for messages.
    std::string logCategory = "ArtHandleTrackerManager";
    
  }; // Config_t
  
  
  /**
   * @brief Constructs a handle manager.
   * @param config (optional) the configuration to be used
   * @see `useEvent()`
   * 
   * An event must be later assigned to it (`useEvent()`) to make it functional.
   */
  ArtHandleTrackerManager(Config_t config = {});
  
  /**
   * @brief Constructs a handle manager covering `event`.
   * @param event the first data viewer (event) object to operate on
   */
  ArtHandleTrackerManager(Event_t const& event, Config_t config = {});
  
  
  // --- BEGIN -- Queries ------------------------------------------------------
  ///@name Queries
  ///@{
  
  /// Returns the number of handles currently tracked.
  unsigned int nTrackedHandles() const;
  
  /// Returns whether the object is associated to an event.
  bool hasEvent() const;
  
  
  /// @}
  // --- END ---- Queries ------------------------------------------------------
  
  
  // --- BEGIN -- Registration of data products --------------------------------
  /// @name Registration of data products
  /// @{
  
  /**
   * @brief Retrieves an handle from the event, and registers it.
   * @tparam T type of data product to retrieve
   * @tparam Args types of the arguments needed by `art::getHandle()`
   * @param args all the arguments that `art::getHandle()` requires
   * @return the handle just read and being managed
   * @see `getValidHandle()`, `registerHandle()`
   * 
   * This function wraps `art::Event::getHandle()`, calling it to obtain the
   * handle and then registering it (like with `registerHandle()`).
   */
  template <typename T, typename... Args>
  art::Handle<T> getHandle(Args&&... args);
  
  /**
   * @brief Retrieves a valid handle from the event, and registers it.
   * @tparam T type of data product to retrieve
   * @tparam Args types of the arguments needed by `art::getValidHandle()`
   * @param args all the arguments that `art::getValidHandle()` requires
   * @return the handle just read and being managed
   * @see `getHandle()`, `registerHandle()`
   * 
   * This is the `art::ValidHandle` sibling of `getHandle()`.
   * See that one for details.
   */
  template <typename T, typename... Args>
  art::ValidHandle<T> getValidHandle(Args&&... args);
  
  /**
   * @brief Registers an existing handle.
   * @tparam Handle the type of handle to register
   * @param handle the handle to be registered
   * @return `handle` (pass through)
   * 
   * This method registers a copy of `handle` into the manager.
   */
  template <typename Handle>
  decltype(auto) registerHandle(Handle&& handle);
  
  /// @}
  // --- END ---- Registration of data products --------------------------------
  
  
  // --- BEGIN -- Operations ---------------------------------------------------
  /// @name Operations
  /// @{
  
  /**
   * @brief Changes the event being served.
   * @param event the new event being served
   * @throw art::Exception (code: `art::errors::LogicError`) if there are still
   *        registered handles
   * 
   * The object starts tracking handles of a new event.
   * 
   * This method requires that any pending handle has been taken care of
   * (even if the new event happens to be the same as the old one).
   * Common options are `removeCachedProductsAndForget()` if cache removal is
   * desired, or `forgetAllHandles()` if it's not.
   */
  void useEvent(Event_t const& event);
  
  
  /**
   * @brief Clears the cached data products for all tracked handles.
   * @return the number of tracked handles which had their cache removed
   * 
   * This method calls `Event_t::removeCachedProduct()` for all tracked handles.
   * This is the core functionality of the manager, which removes the cache
   * of all tracked handles.
   * 
   * The _art_ framework always makes the handles used to remove the cache
   * invalid (`Handle::clear()`).
   * After the removal, the object stops tracking the handles (like with a call
   * to `forgetAllHandles()`).
   * 
   * Calling this method when there are no handles registered has no effect.
   */
  unsigned int removeCachedProducts();
  
  /// Stops tracking any handle, without removing their cache first.
  void forgetAllHandles();
  
  
  /**
   * @brief Completes the work on the associated `event`.
   * @param removeCache (default: `true`) removes tracked data product caches
   * @param event if specified, only acts on cache from this `event`
   * @return the number of cleared data product caches
   * 
   * This is a shortcut for having optional release of cache in a single call,
   * depending on the value of `removeCache`:
   *  * if `true`, `removeCachedProducts()` is called and its count is returned;
   *  * if `false`, `forgetAllHandles()` is called and `0` is returned.
   * 
   * In addition, the object is disassociated from the event, and a call to
   * `useEvent()` will be needed before this object can be operative again.
   * 
   * If the `event` parameter is not null, a check is done that the event being
   * currently used matches `event`, and if not no operation happens (`0` is
   * returned, but the object is not disassociated from its current event).
   * 
   */
  unsigned int doneWithEvent
    (bool removeCache = true, art::Event const* event = nullptr);
  
  /// @}
  // --- END ---- Operations ---------------------------------------------------
  
  
    private:
  
  /// Type of pointer to any handle tracker.
  using TrackerPtr = std::unique_ptr<util::ArtHandleTrackerInterface<Event_t>>;
  
  
  // --- BEGIN -- Data ---------------------------------------------------------
  Event_t const* fEvent = nullptr; ///< Event being manager. Must be present!
  
  Config_t const fConfig; ///< Configuration.
  
  
  std::vector<TrackerPtr> fTrackers; ///< List of managed handle trackers.
  
  // --- END ---- Data ---------------------------------------------------------
  
  
  /**
   * @brief Returns the pointer to an handle equivalent to `like`.
   * @tparam Handle the type of the handle being queried
   * @return a pointer to the handle, or `nullptr` if not found
   * 
   * This is an half-assed method since it needs to know already the exact type
   * of the handle being queried (it's not even agnostic to the difference
   * between `art::Handle` and `art::ValidHandle`).
   * For this reason, it's not really useful to the users, who would probably
   * want to know if the data product is cached, via any mean.
   * 
   * 
   * ### Technical note
   * 
   * Delegating the matching to `util::ArtHandleTrackerInterface` is not possible
   * (`like` is a template parameter and can't be passed via virtual interface)
   * and even working around that and finding the match, the returned value needs
   * to be of type `Handle`, which is not necessarily the type of handle stored
   * in the tracker.
   * 
   * For the full functionality of knowing if a data product is tracked, a
   * separate registry may be kept, still complicated by the fact that part of
   * the registry entry is a C++ type (we may need to store a sanitized
   * `std::type_info` for that).
   */
  template <typename Handle>
  Handle const* findHandle(Handle const& like) const;
  
  
  /// Returns whether `tracker` tracks the same product that `handle` handles.
  template <typename Handle>
  static bool handleSameProduct(TrackerPtr tracker, Handle const& handle);
  
  
  /// Checks that the object is in a state where it can perform operations.
  /// @param where identification of the calling function for error messages
  /// @throw art::Exception (code: `art::errors::LogicError`) on failure
  void canOperate(const char* where) const;
  
  
}; // util::ArtHandleTrackerManager


// -----------------------------------------------------------------------------
/**
 * @brief Variant of `util::ArtHandleTrackerManager` in local scope.
 * @tparam Event the type of data product source (the "principal") to serve
 *
 * This class allows the removal of data products cached in an event.
 * The general functionality is described in `util::ArtHandleTrackerManager`.
 * 
 * This class in particular offers a tuned interface for the case where the
 * object is local (which is convenient for multithreading and if all the
 * reading and usage happens in the same scope):
 *  * the object supports being associated to only one event in its lifetime,
 *    and the association must be established on construction;
 *  * removal is automatically performed on destruction at the exit of the
 *    scope of the object (this is "Resource Acquisition Is Initialization"
 *    idiom)
 * 
 * Whether the object will actually remove the caches can be decided on
 * construction (by default it does), and then changed with
 * `setRemoveCachedProducts()`.
 * The removal can also be anticipated by explicitly calling `doneWithEvent()`
 * (no `event` parameter is supported), after which the destruction will not
 * do anything.
 *
 * Example:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * void MyModule::analyze(art::Event const & event) override {
 * 
 *   util::LocalArtHandleTrackerManager dataCacheRemover{ event };
 *   
 *   auto const& handle = dataCacheRemover.getHandle<MyDataProduct>(fTag);
 *   
 *   auto results = processData(*handle); // ... or whatever
 *   
 *   results.Write(); // ... or another portion of whatever
 * }
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 */
template <typename Event>
class util::LocalArtHandleTrackerManager {
  
  /// The actual manager doing the work.
  ArtHandleTrackerManager<Event> fManager;

  bool fRemoveCache; ///< Whether to remove cache on destructor.
  
    public:
  
  /// Type of the _art_ event.
  using Event_t = typename ArtHandleTrackerManager<Event>::Event_t;
  
  /**
   * @brief Constructor: operates on the specified `event`.
   * @param event the event this object will operate on
   * @param removeCache (default: `true`) remove tracked data caches when done
   * 
   * 
   * The parameter `removeCache` controls whether the data caches are removed
   * on destruction. It may be changed along the way with
   * `setRemoveCachedProducts()`.
   * 
   */
  LocalArtHandleTrackerManager(Event const& event, bool removeCache = true);
  
  /// Destructor; will implement the RAII pattern (i.e. call `doneWithEvent()`).
  ~LocalArtHandleTrackerManager();
  
  
  // --- BEGIN -- Queries ------------------------------------------------------
  ///@name Queries
  ///@{
  
  /// Returns the number of handles currently tracked.
  unsigned int nTrackedHandles() const;
  
  /// Returns whether caches will be removed on destruction.
  /// @see `setRemoveCachedProducts()`
  bool willRemoveCachedProducts() const;
  
  /// @}
  // --- END ---- Queries ------------------------------------------------------
  
  
  // --- BEGIN -- Registration of data products --------------------------------
  /// @name Registration of data products
  /// @{
  
  /**
   * @brief Retrieves an handle from the event, and registers it.
   * @tparam T type of data product to retrieve
   * @tparam Args types of the arguments needed by `art::getHandle()`
   * @param args all the arguments that `art::getHandle()` requires
   * @return the handle just read and being managed
   * @see `getValidHandle()`, `registerHandle()`
   * 
   * This function wraps `art::Event::getHandle()`, calling it to obtain the
   * handle and then registering it (like with `registerHandle()`).
   */
  template <typename T, typename... Args>
  art::Handle<T> getHandle(Args&&... args);
  
  /**
   * @brief Retrieves a valid handle from the event, and registers it.
   * @tparam T type of data product to retrieve
   * @tparam Args types of the arguments needed by `art::getValidHandle()`
   * @param args all the arguments that `art::getValidHandle()` requires
   * @return the handle just read and being managed
   * @see `getHandle()`, `registerHandle()`
   * 
   * This is the `art::ValidHandle` sibling of `getHandle()`.
   * See that one for details.
   */
  template <typename T, typename... Args>
  art::ValidHandle<T> getValidHandle(Args&&... args);
  
  /**
   * @brief Registers an existing handle.
   * @tparam Handle the type of handle to register
   * @param handle the handle to be registered
   * @return the same `handle` is passed through
   * 
   * This method registers a copy of `handle` into the manager.
   */
  template <typename Handle>
  decltype(auto) registerHandle(Handle&& handle);
  
  /// @}
  // --- END ---- Registration of data products --------------------------------
  
  
  // --- BEGIN -- Operations ---------------------------------------------------
  /// @name Operations
  /// @{
  
  /// @brief Sets whether to remove tracked data caches on destruction or not.
  /// @param removeCache if `true`, data caches will be removed on destruction
  void setRemoveCachedProducts(bool removeCache);
  
  /**
   * @brief Clears the cached data products for all tracked handles.
   * @return the number of tracked handles which had their cache removed
   * 
   * This method calls `Event_t::removeCachedProduct()` for all tracked handles.
   * This is the core functionality of the manager, which removes the cache
   * of all tracked handles.
   * 
   * The _art_ framework always makes the handles used to remove the cache
   * invalid (`Handle::clear()`).
   * After the removal, the object stops tracking the handles (like with a call
   * to `forgetAllHandles()`).
   * 
   * Calling this method when there are no handles registered has no effect.
   */
  unsigned int removeCachedProducts();
  
  /// Stops tracking any handle, without removing their cache first.
  void forgetAllHandles();
  
  
  /**
   * @brief Completes the work on the associated `event`.
   * @param removeCache (default: `true`) removes tracked data product caches
   * @return the number of cleared data product caches
   * @see `doneWithEvent()`
   * 
   * This is a shortcut for having optional release of cache in a single call,
   * depending on the value of `removeCache`:
   *  * if `true`, `removeCachedProducts()` is called and its count is returned;
   *  * if `false`, `forgetAllHandles()` is called and `0` is returned.
   * 
   * Differently from `util::ArtHandleTrackerManager`, the object is _not_
   * disassociated from the event: new handles can be registered and the
   * clearing settings are preserved as they are.
   * To be really done with the event, `setRemoveCachedProducts()` needs to be
   * set to `false`, and users should refrain from calling
   * `removeCachedProducts().
   * 
   * This method ignores and overrides the value set by
   * `setRemoveCachedProducts()` and the one at construction.
   */
  unsigned int doneWithEvent(bool removeCache);
  
  
  /**
   * @brief Completes the work on the associated `event`.
   * @return the number of cleared data product caches
   * @see `doneWithEvent(bool)`, `setRemoveCachedProducts()`
   * 
   * This is a shortcut for having optional release of cache in a single call,
   * depending on the value of `willRemoveCachedProducts()`:
   *  * if `true`, `removeCachedProducts()` is called and its count is returned;
   *  * if `false`, `forgetAllHandles()` is called and `0` is returned.
   * 
   * In addition, the object is disassociated from the event, and the object
   * will not be available for use any more.
   */
  unsigned int doneWithEvent();
  
  /// @}
  // --- END ---- Operations ---------------------------------------------------
  
}; // util::LocalArtHandleTrackerManager


// -----------------------------------------------------------------------------
/**
 * @brief Interface to facilitate the use of `util::ArtHandleTracker`
 *        specializations.
 * 
 * This is NOT able to return the type of handle it's handling.
 * 
 */
template <typename Event>
struct util::ArtHandleTrackerInterface {
  
  /// Virtual destructor (does nothing).
  virtual ~ArtHandleTrackerInterface() = default;
  
  /// Removes the cached data product from `event`.
  /// Handle is cleared and won't be valid any more.
  bool removeCachedProduct() { return doRemoveCachedProduct(); }
  
  
  /// Returns a container for a pointer to the managed handle. Caveat emptor.
  std::any handlePtr() const { return doHandlePtr(); }
  
  /// Returns the name of the class of handled data product.
  std::string productClass() const { return doProductClass(); }
  
  /// Returns the name of the class of handled data product.
  std::type_info const* productType() const { return doProductType(); }
  
  /// Returns the tag of handled data product.
  art::InputTag inputTag() const { return doInputTag(); }
  
  /// Returns whether this and the `other` objects handle the same data product.
  bool hasSameDataProduct
    (util::ArtHandleTrackerInterface<Event> const& other) const
    {
      return inputTag() == other.inputTag()
        && productType() == other.productType();
    }
  
    protected:
  
  /// Deferred implementation of `removeCachedProduct()`.
  virtual bool doRemoveCachedProduct() = 0;
  
  /// Deferred implementation of `handlePtr()`.
  virtual std::any doHandlePtr() const = 0;
  
  /// Deferred implementation of `productClass()`.
  virtual std::string doProductClass() const = 0;
  
  /// Deferred implementation of `productType()`.
  virtual std::type_info const* doProductType() const = 0;
  
  /// Deferred implementation of `inputTag()`.
  virtual art::InputTag doInputTag() const = 0;
  
}; // ArtHandleTrackerInterface<>


// -----------------------------------------------------------------------------
// ---  Template implementation
// -----------------------------------------------------------------------------
// ---  util::details::ArtHandleTracker
// -----------------------------------------------------------------------------
namespace util::details {
  template <typename Handle, typename Enable = void> class ProvenanceGetter;
  template <typename T, typename Event> class ArtHandleTracker;
  
  /// use candy
  template <typename Handle>
  std::string productClassOf(Handle const& handle);
  template <typename Handle>
  std::type_info const* productTypeOf(Handle const& handle);
  template <typename Handle>
  art::InputTag inputTagOf(Handle const& handle);
  
} // namespace util::details

// -----------------------------------------------------------------------------
/// Helper to extract basic information from one handle.
/// The default implementation supports `art::Handle` and `art::ValidHandle`.
template <typename Handle, typename>
class util::details::ProvenanceGetter {
  
  static art::Provenance const* provenance(Handle const& handle)
    { return handle.provenance(); }
  
    public:
  /// Returns the name of the class pointed by the handle.
  static std::string productClass(Handle const& handle)
    {
      auto const* prov = provenance(handle);
      return prov? prov->producedClassName(): "";
    }
  
  /// Returns the C++ type of the handled data product.
  static std::type_info const* productType()
    { return &typeid(typename Handle::element_type); }
  
  /// Returns the C++ type of the handled data product.
  static std::type_info const* productType(Handle const&)
    { return productType(); }
  
  /// Returns the input tag of the handled data product.
  /// Deferred implementation of `inputTag()`.
  static art::InputTag inputTag(Handle const& handle)
    {
      auto const* prov = provenance(handle);
      return prov? prov->inputTag(): art::InputTag{};
    }
  
  
}; // util::details::ProvenanceGetter


template <typename Handle>
std::string util::details::productClassOf(Handle const& handle)
  { return ProvenanceGetter<Handle>::productClass(handle); }
  
template <typename Handle>
std::type_info const* util::details::productTypeOf(Handle const& handle)
  { return ProvenanceGetter<Handle>::productType(handle); }
  
template <typename Handle>
art::InputTag util::details::inputTagOf(Handle const& handle)
  { return ProvenanceGetter<Handle>::inputTag(handle); }


// -----------------------------------------------------------------------------
/**
 * @brief Tracks _art_ handle objects.
 * @tparam Handle type of handle being tracked
 * 
 * @note Due to my limited expertise with metaprogramming (and/or for C++
 *       limitations?) I need one tracker per data type.
 *       The `util::ArtHandleTrackerManager` class should mitigate the hassle
 *       deriving from this limitation.
 * 
 */
template <typename Event, typename Handle>
class util::details::ArtHandleTracker: public ArtHandleTrackerInterface<Event> {
  
  Event const* fEvent = nullptr;
  
  Handle fHandle;
  
  // --- BEGIN -- Virtual method implementation --------------------------------
  /// @name Virtual method implementation
  /// @{
  
  /// Actually removes the cache.
  virtual bool doRemoveCachedProduct() override
    {
      mf::LogDebug{ "ArtHandleTracker" } << "Removing cache for handle<"
        << fHandle.provenance()->producedClassName()
        << ">(" << fHandle.provenance()->inputTag().encode() << ").";
      return fEvent->removeCachedProduct(fHandle);
    }
  
  /// Returns a pointer to the managed handle, wrapped in `std::any`.
  virtual std::any doHandlePtr() const override
    { return { &handle() }; }
  
  /// Returns the name of the class pointed by the handle.
  virtual std::string doProductClass() const override
    { return productClassOf(fHandle); }
  
  /// Returns the C++ type of the handled data product.
  virtual std::type_info const* doProductType() const override
    { return ProvenanceGetter<Handle>::productType(); }
  
  /// Returns the input tag of the handled data product.
  /// Deferred implementation of `inputTag()`.
  virtual art::InputTag doInputTag() const override
    { return inputTagOf(fHandle); }

  /// @}
  // --- END ---- Virtual method implementation --------------------------------
  
    public:
  
  /// Constructor: records all the needed information.
  ArtHandleTracker(Event const& event, Handle handle)
    : fEvent(&event), fHandle(std::move(handle)) {}
  
  /// Returns the tracked handle.
  Handle const& handle() const { return fHandle; }
  
  /// Returns the provenance information of the handle.
  art::Provenance const* provenance() const { return fHandle.provenance(); }
  
}; // util::details::ArtHandleTracker<>


// -----------------------------------------------------------------------------
// --- util::ArtHandleTrackerManager
// -----------------------------------------------------------------------------
template <typename Event>
util::ArtHandleTrackerManager<Event>::ArtHandleTrackerManager
  (Config_t config /* = {} */)
  : fConfig{ std::move(config) }
{}


// -----------------------------------------------------------------------------
template <typename Event>
util::ArtHandleTrackerManager<Event>::ArtHandleTrackerManager
  (Event_t const& event, Config_t config /* = {} */)
  : fEvent{ &event }, fConfig{ std::move(config) }
{}


// -----------------------------------------------------------------------------
template <typename Event>
void util::ArtHandleTrackerManager<Event>::useEvent(Event_t const& event) {
  
  if (nTrackedHandles() > 0) {
    // since fEvent might be invalid, we don't attempt to figure out which ID
    // that event might have had
    throw art::Exception{ art::errors::LogicError }
      << "ArtHandleTrackerManager attempted to change event to "
      << event.id() << " when " << nTrackedHandles()
      << " handles are still tracked.\n";
  }
  
  fEvent = &event;
  
} // util::ArtHandleTrackerManager<Event>::useEvent()


// -----------------------------------------------------------------------------
template <typename Event>
template <typename T, typename... Args>
art::Handle<T> util::ArtHandleTrackerManager<Event>::getHandle
  (Args&&... args)
{
  canOperate("getHandle()");
  return registerHandle
    (fEvent->template getHandle<T>(std::forward<Args>(args)...));
}


// -----------------------------------------------------------------------------
template <typename Event>
template <typename T, typename... Args>
art::ValidHandle<T> util::ArtHandleTrackerManager<Event>::getValidHandle
  (Args&&... args)
{
  canOperate("getValidHandle()");
  return registerHandle
    (fEvent->template getValidHandle<T>(std::forward<Args>(args)...));
}


// -----------------------------------------------------------------------------
template <typename Event>
template <typename Handle>
auto util::ArtHandleTrackerManager<Event>::registerHandle
  (Handle&& handle) -> decltype(auto)
{
  using util::details::ProvenanceGetter;
  
  using Handle_t = std::decay_t<Handle>;
  
  canOperate("registerHandle()");
  
  using Tracker_t = util::details::ArtHandleTracker<Event_t, Handle_t>;
  
  // if it's already registered, we don't want to have it again
  if (auto ptr = findHandle(handle)) {
    auto const& registeredHandle = *ptr;
    mf::LogDebug msg { fConfig.logCategory };
    msg
      << "Handle<" << details::productClassOf(registeredHandle)
      << ">(" << details::inputTagOf(registeredHandle).encode()
      << ") was already registered";
    if (typeid(registeredHandle) != typeid(Handle_t)) {
      msg << " as a different handle type (" << typeid(registeredHandle).name()
        << ", now " << typeid(Handle_t).name() << ")";
    }
    msg << ".";
  }
  else {
    mf::LogDebug{ fConfig.logCategory }
      << "Registering handle<" << details::productClassOf(handle)
      << ">(" << details::inputTagOf(handle).encode()
      << ") (handle type: " << typeid(Handle_t).name() << ")."
      ;
    
    fTrackers.push_back(std::make_unique<Tracker_t>(*fEvent, handle));
  }
  
  return std::forward<Handle>(handle);
  
} // util::ArtHandleTrackerManager<>::registerHandle()


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::ArtHandleTrackerManager<Event>::nTrackedHandles() const
  { return fTrackers.size(); }


// -----------------------------------------------------------------------------
template <typename Event>
bool util::ArtHandleTrackerManager<Event>::hasEvent() const
  { return fEvent; }


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::ArtHandleTrackerManager<Event>::removeCachedProducts() {
  
  // we remove cache in opposite order to registration: it's a C++ tradition.
  unsigned int const nRemoved = std::count_if(
    fTrackers.crbegin(), fTrackers.crend(),
    [](auto const& tracker){ return tracker->removeCachedProduct(); }
    );
  
  forgetAllHandles();
  
  return nRemoved;
} // util::ArtHandleTrackerManager<Event>::removeCachedProducts()


// -----------------------------------------------------------------------------
template <typename Event>
void util::ArtHandleTrackerManager<Event>::forgetAllHandles()
  { fTrackers.clear(); }

  


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::ArtHandleTrackerManager<Event>::doneWithEvent
  (bool removeCache /* = true */, art::Event const* event /* = nullptr */)
{
  if (event && (fEvent != event)) return 0;
  
  unsigned int count = 0;
  if (removeCache) count = removeCachedProducts();
  else count = forgetAllHandles();
  
  fEvent = nullptr;
  
  return count;
} // util::ArtHandleTrackerManager<Event>::doneWithEvent()


// -----------------------------------------------------------------------------
template <typename Event>
template <typename Handle>
Handle const* util::ArtHandleTrackerManager<Event>::findHandle
  (Handle const& like) const
{
  
  // look for it, one by one
  for (TrackerPtr const& tracker: fTrackers) {
    
    
    
    std::any anyHandlePtr = tracker->handlePtr();
    Handle const** handlePtr = std::any_cast<Handle const*>(&anyHandlePtr);
    if (!handlePtr) continue; // not the right type
    
    if ((*handlePtr)->provenance()->inputTag() != like.provenance()->inputTag())
      continue; // different tag
    
    return *handlePtr;
  } // for
  
  return nullptr; // failed
  
} // util::ArtHandleTrackerManager<Event>::findTracker()


// -----------------------------------------------------------------------------
template <typename Event>
template <typename Handle>
bool util::ArtHandleTrackerManager<Event>::handleSameProduct
  (TrackerPtr tracker, Handle const& handle)
{
  using util::details::ProvenanceGetter;
  
  if (!tracker) return false;
  if (ProvenanceGetter<Handle>::productType() != tracker->productType())
    return false;
  if (ProvenanceGetter<Handle>::inputTag() != tracker->inputTag())
    return false;
  
  return true;
} // util::ArtHandleTrackerManager<>::handleSameProduct()


// -----------------------------------------------------------------------------
template <typename Event>
void util::ArtHandleTrackerManager<Event>::canOperate(const char* where) const {
  
  if (fEvent) return;
  
  throw art::Exception(art::errors::LogicError)
    << "util::ArtHandleTrackerManager attempted an operation without a valid"
      " event.\nThe operation was: " << where << "\n";
  
} // util::ArtHandleTrackerManager<Event>::canOperate()


// -----------------------------------------------------------------------------
// ---  util::LocalArtHandleTrackerManager
// -----------------------------------------------------------------------------
template <typename Event>
util::LocalArtHandleTrackerManager<Event>::LocalArtHandleTrackerManager
  (Event const& event, bool removeCache /* = true */)
  : fManager{ event }, fRemoveCache{ removeCache }
  {}


// -----------------------------------------------------------------------------
template <typename Event>
util::LocalArtHandleTrackerManager<Event>::~LocalArtHandleTrackerManager()
  { doneWithEvent(); }


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::LocalArtHandleTrackerManager<Event>::nTrackedHandles() const
  {return fManager.nTrackedHandles(); }


// -----------------------------------------------------------------------------
template <typename Event>
bool util::LocalArtHandleTrackerManager<Event>::willRemoveCachedProducts() const
  { return fRemoveCache; }


// -----------------------------------------------------------------------------
template <typename Event>
template <typename T, typename... Args>
art::Handle<T> util::LocalArtHandleTrackerManager<Event>::getHandle
  (Args&&... args)
  { return fManager.template getHandle<T>(std::forward<Args>(args)...); }


// -----------------------------------------------------------------------------
template <typename Event>
template <typename T, typename... Args>
art::ValidHandle<T> util::LocalArtHandleTrackerManager<Event>::getValidHandle
  (Args&&... args)
  { return fManager.template getValidHandle<T>(std::forward<Args>(args)...); }


// -----------------------------------------------------------------------------
template <typename Event>
template <typename Handle>
auto util::LocalArtHandleTrackerManager<Event>::registerHandle
  (Handle&& handle) -> decltype(auto)
{
  return fManager.template registerHandle<Handle>(std::forward<Handle>(handle)); 
}


// -----------------------------------------------------------------------------
template <typename Event>
void util::LocalArtHandleTrackerManager<Event>::setRemoveCachedProducts
  (bool removeCache)
  { fRemoveCache = removeCache; }


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::LocalArtHandleTrackerManager<Event>::removeCachedProducts()
  { return fManager.removeCachedProducts(); }


// -----------------------------------------------------------------------------
template <typename Event>
void util::LocalArtHandleTrackerManager<Event>::forgetAllHandles()
  { return fManager.forgetAllHandles(); }


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::LocalArtHandleTrackerManager<Event>::doneWithEvent
  (bool removeCache)
{
  assert(fManager.hasEvent());
  
  if (removeCache) return removeCachedProducts();
  
  forgetAllHandles();
  return 0;
} // util::LocalArtHandleTrackerManager<Event>::doneWithEvent()


// -----------------------------------------------------------------------------
template <typename Event>
unsigned int util::LocalArtHandleTrackerManager<Event>::doneWithEvent()
  { return doneWithEvent(willRemoveCachedProducts()); }


// -----------------------------------------------------------------------------


#endif // SBNDCORE_DECODERS_DECODERTOOLS_DUMPERS_ARTHANDLETRACKERMANAGER_H
