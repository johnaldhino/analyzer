//*-- Author :    Ole Hansen   16/05/00

//////////////////////////////////////////////////////////////////////////
//
// THaDetMap
//
// The standard detector map for a Hall A detector.
//
//////////////////////////////////////////////////////////////////////////

#include "THaDetMap.h"
#include "Decoder.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Decoder;

#define ALL(c) (c).begin(), (c).end()

// "Database" of known frontend module numbers and types
// FIXME: load from db_cratemap
struct ModuleDef {
  UInt_t      model; // model identifier
  ChannelType type;  // Module type
} __attribute__((aligned(8)));

static const vector<ModuleDef> module_list {
  { 1875, ChannelType::kTDC },
  { 1877, ChannelType::kTDC },
  { 1881, ChannelType::kADC },
  { 1872, ChannelType::kTDC },
  { 3123, ChannelType::kADC },
  { 1182, ChannelType::kADC },
  { 792,  ChannelType::kADC },
  { 775,  ChannelType::kTDC },
  { 767,  ChannelType::kTDC },
  { 3201, ChannelType::kTDC },
  { 6401, ChannelType::kTDC },
  { 1190, ChannelType::kTDC },
  { 250,  ChannelType::kMultiFunctionADC },
};

//_____________________________________________________________________________
void THaDetMap::Module::SetResolution( Double_t res )
{
  resolution = res;
}

//_____________________________________________________________________________
void THaDetMap::Module::SetTDCMode( Bool_t cstart )
{
  cmnstart = cstart;
  if( cstart && type == ChannelType::kCommonStopTDC )
    type = ChannelType::kTDC;
  else if( !cstart && type == ChannelType::kTDC )
    type = ChannelType::kCommonStopTDC;
}

//_____________________________________________________________________________
void THaDetMap::Module::MakeTDC()
{
  if( IsCommonStart() )
    type = ChannelType::kTDC;
  else
    type = ChannelType::kCommonStopTDC;
}

//_____________________________________________________________________________
void THaDetMap::Module::MakeADC()
{
  type = ChannelType::kADC;
}

//_____________________________________________________________________________
void THaDetMap::CopyMap( const ModuleVec_t& map )
{
  // Deep-copy the vector of module pointers
  fMap.reserve(map.capacity());
  for( const auto& m : map ) {
    fMap.emplace_back(
#if __cplusplus >= 201402L
      std::make_unique<Module>(*m)
#else
      new Module(*m)
#endif
    );
  }
}

//_____________________________________________________________________________
THaDetMap::THaDetMap( const THaDetMap& rhs )
{
  // Copy constructor

  CopyMap(rhs.fMap);
}

//_____________________________________________________________________________
THaDetMap& THaDetMap::operator=( const THaDetMap& rhs )
{
  // THaDetMap assignment operator

  if( this != &rhs ) {
    fMap.clear();
    CopyMap(rhs.fMap);
  }
  return *this;
}

//_____________________________________________________________________________
Int_t THaDetMap::AddModule( UShort_t crate, UShort_t slot,
                            UShort_t chan_lo, UShort_t chan_hi,
                            UInt_t first, UInt_t model, Int_t refindex,
                            Int_t refchan, UInt_t plane, UInt_t signal )
{
  // Add a module to the map.

  // Logical channels can run either "with" or "against" the physical channel
  // numbers:
  // lo<=hi :    lo -> first
  //             hi -> first+hi-lo
  //
  // lo>hi  :    hi -> first
  //             lo -> first+lo-hi
  //
  // To indicate the second case, The flag "reverse" is set to true in the
  // module. The fields lo and hi are reversed so that lo<=hi always.
  //
  bool reverse = (chan_lo > chan_hi);

#if __cplusplus >= 201402L
  auto pm = make_unique<Module>();
#else
  auto pm = unique_ptr<Module>(new THaDetMap::Module);
#endif
  Module& m = *pm;
  m.crate = crate;
  m.slot  = slot;
  if( reverse ) {
    m.lo = chan_hi;
    m.hi = chan_lo;
  } else {
    m.lo = chan_lo;
    m.hi = chan_hi;
  }
  m.first = first;
  m.model = model;
  m.refindex = refindex;
  m.refchan  = refchan;
  m.plane = plane;
  m.signal = signal;
  m.SetResolution(0.0);
  m.reverse = reverse;
  m.cmnstart = false;

  // Set the module type using our lookup table of known model numbers.
  // If the model is not predefined, this can be done manually later with
  // calls to MakeADC()/MakeTDC().
  if( model != 0 ) {
    auto it = find_if(ALL(module_list),
                      [model]( const ModuleDef& m ) { return m.model == model; });
    m.type = (it != module_list.end()) ? it->type : ChannelType::kUndefined;
    if( m.type == ChannelType::kTDC )
      m.type = ChannelType::kCommonStopTDC; // common stop is the default
  } else
    m.type = ChannelType::kUndefined;

  fMap.push_back(std::move(pm));
  return fMap.size();
}

//_____________________________________________________________________________
THaDetMap::Module* THaDetMap::Find( UShort_t crate, UShort_t slot,
                                    UShort_t chan )
{
  // Return the module containing crate, slot, and channel chan.
  // If several matching modules exist (which mean the map is misconfigured),
  // only the first one is returned. If none match, return nullptr.
  // Since the map is usually small and not necessarily sorted, a simple
  // linear search is done.

  auto found = find_if( ALL(fMap), [crate,slot,chan](const unique_ptr<Module>& d)
  { return ( d->crate == crate && d->slot == slot &&
             d->lo <= chan && chan <= d->hi );
  });

  return (found != fMap.end()) ? found->get() : nullptr;
}

//_____________________________________________________________________________
Int_t THaDetMap::Fill( const vector<Int_t>& values, UInt_t flags )
{
  // Fill the map with 'values'. Depending on 'flags', the values vector
  // is interpreted as a 4-, 5-, 6- or 7-tuple:
  //
  // The first 4 values are interpreted as (crate,slot,start_chan,end_chan)
  // Each of the following flags causes one more value to be used as part of
  // the tuple for each module:
  //
  // kFillLogicalChannel - Logical channel number for 'start_chan'.
  // kFillModel          - The module's hardware model number (see AddModule())
  // kFillRefChan        - Reference channel (for pipeline TDCs etc.)
  // kFillRefIndex       - Reference index (for pipeline TDCs etc.)
  // kFillPlane          - Which plane in detector (for Hall C)
  // kFillSignal         - Which signal type (for Hall C)
  //
  // If more than one flag is present, the numbers will be interpreted
  // in the order the flags are listed above.
  // Example:
  //      flags = kFillModel | kFillRefChan
  // ==>
  //      the vector is interpreted as a series of 6-tuples in the order
  //      (crate,slot,start_chan,end_chan,model,refchan).
  //
  // If kFillLogicalChannel is not set then the first logical channel numbers
  // are automatically calculated for each module, assuming the numbers are
  // sequential.
  //
  // By default, an existing map is overwritten. If the flag kDoNotClear 
  // is present, then the data are appended.
  //
  // The return value is the number of modules successfully added,
  // or negative if an error occurred.

  typedef vector<Int_t>::size_type vsiz_t;

  if( (flags & kDoNotClear) == 0 )
    Clear();

  vsiz_t tuple_size = 4;
  if( flags & kFillLogicalChannel )
    tuple_size++;
  if( flags & kFillModel )
    tuple_size++;
  if( flags & kFillRefChan )
    tuple_size++;
  if( flags & kFillRefIndex )
    tuple_size++;
  if( flags & kFillPlane )
    tuple_size++;
  if( flags & kFillSignal )
    tuple_size++;

  // For historical reasons, logical channel numbers start at 1
  UInt_t prev_first = 1, prev_nchan = 0;
  // Defaults for optional values
  UInt_t first = 0, model = 0;
  Int_t rchan = -1, ref = -1;
  UInt_t plane = 0, signal = 0;

  Int_t ret = 0;
  for( vsiz_t i = 0; i < values.size(); i += tuple_size ) {
    // For compatibility with older maps, crate < 0 means end of data
    if( values[i] < 0 )
      break;
    // Now we require a full tuple
    if( i + tuple_size > values.size() ) {
      ret = -127;
      break;
    }

    vsiz_t k = 4;
    if( flags & kFillLogicalChannel )
      first = values[i + k++];
    else {
      first = prev_first + prev_nchan;
    }
    if( flags & kFillModel )
      model = values[i + k++];
    if( flags & kFillRefChan )
      rchan = values[i + k++];
    if( flags & kFillRefIndex )
      ref = values[i + k++];
    if( flags & kFillPlane )
      plane = values[i + k++];
    if( flags & kFillSignal )
      signal = values[i + k++];

    ret = AddModule(values[i], values[i + 1], values[i + 2], values[i + 3],
                    first, model, ref, rchan, plane, signal);
    if( ret <= 0 )
      break;
    prev_first = first;
    prev_nchan = GetNchan(ret - 1);
  }

  return ret;
}

//_____________________________________________________________________________
Int_t THaDetMap::GetTotNumChan() const
{
  // Get sum of the number of channels of all modules in the map. This is
  // typically the total number of hardware channels used by the detector.

  Int_t sum = 0;
  for( const auto & m : fMap )
    sum += m->GetNchan();

  return sum;
}


//_____________________________________________________________________________
void THaDetMap::GetMinMaxChan( Int_t& min, Int_t& max, ECountMode mode ) const
{
  // Put the minimum and maximum logical or reference channel numbers
  // into min and max. If refidx is true, check refindex, else check logical
  // channel numbers.

  min = kMaxInt;
  max = kMinInt;
  bool do_ref = (mode == kRefIndex);
  for( const auto& m : fMap ) {
    Int_t m_min = do_ref ? m->refindex : m->first;
    Int_t m_max = do_ref ? m->refindex : m->first + m->hi - m->lo;
    if( m_min < min )
      min = m_min;
    if( m_max > max )
      max = m_max;
  }
}


//_____________________________________________________________________________
void THaDetMap::Print( Option_t* ) const
{
  // Print the contents of the map

  cout << "Size: " << GetSize() << endl;
  for( const auto& m : fMap ) {
    cout << " "
         << setw(5) << m->crate
         << setw(5) << m->slot
         << setw(5) << m->lo
         << setw(5) << m->hi
         << setw(5) << m->first
         << setw(5) << m->model;
    if( m->IsADC() )
      cout << setw(4) << " ADC";
    if( m->IsTDC() )
      cout << setw(4) << " TDC";
    cout << setw(5) << m->refchan
         << setw(5) << m->refindex
         << setw(8) << m->resolution
         << setw(5) << m->plane
         << setw(5) << m->signal
         << endl;
  }
}

//_____________________________________________________________________________
void THaDetMap::Reset()
{
  // Clear() the map and reset the array size to zero, freeing memory.

  Clear();
  fMap.shrink_to_fit(); // FWIW
}

//_____________________________________________________________________________
void THaDetMap::Sort()
{
  // Sort the map by crate/slot/low channel

  sort( ALL(fMap),
	[]( const unique_ptr<Module>& a, const unique_ptr<Module>& b) {
    if( a->crate < b->crate ) return true;
    if( a->crate > b->crate ) return false;
    if( a->slot  < b->slot )  return true;
    if( a->slot  > b->slot )  return false;
    return (a->lo < b->lo);
  });

}

//=============================================================================
THaDetMap::Iterator
THaDetMap::MakeIterator( const THaEvData& evdata )
{
  return THaDetMap::Iterator(*this, evdata);
}

//_____________________________________________________________________________
THaDetMap::MultiHitIterator
THaDetMap::MakeMultiHitIterator( const THaEvData& evdata )
{
  return THaDetMap::MultiHitIterator(*this, evdata);
}

//_____________________________________________________________________________
THaDetMap::Iterator::Iterator( THaDetMap& detmap, const THaEvData& evdata )
  : fDetMap(detmap), fEvData(evdata), fMod(nullptr), fNMod(fDetMap.GetSize()),
    fNTotChan(fDetMap.GetTotNumChan()), fNChan(0), fIMod(-1), fIChan(-1)
{
  fHitInfo.ev = fEvData.GetEvNum();
  // Initialize iterator state to point to the first item
  ++(*this);
}

//_____________________________________________________________________________
string THaDetMap::Iterator::msg( const char* txt ) const
{
  // Format message string for exceptions
  ostringstream ostr;
  ostr << "Event " << fHitInfo.ev << ", ";
  ostr << "channel "
       << fHitInfo.crate << "/" << fHitInfo.slot << "/"
       << fHitInfo.chan  << "/" << fHitInfo.hit << ": ";
  if( txt and *txt ) ostr << txt; else ostr << "Unspecified error.";
  return ostr.str();
}

#define CRATE_SLOT( x ) (x).crate, (x).slot
static inline bool LOOPDONE( Int_t i, Int_t n ) { return i < 0 or i >= n; }
static inline bool NO_NEXT( Int_t& i, Int_t n ) { return i >= n or ++i >= n; }

//_____________________________________________________________________________
THaDetMap::Iterator& THaDetMap::Iterator::operator++()
{
  // Advance iterator to next active channel
 nextmod:
  if( LOOPDONE(fIChan, fNChan) ) {
    if( NO_NEXT(fIMod, fNMod) ) {
      fHitInfo.reset();
      return *this;
    }
    fMod = fDetMap.GetModule(fIMod);
    if( !fMod )
      throw std::logic_error("NULL detector map module. Program bug. Call expert.");
    fHitInfo.set_crate_slot(fMod);
    fHitInfo.module = fEvData.GetModule(CRATE_SLOT(fHitInfo));
    fNChan = fEvData.GetNumChan(CRATE_SLOT(fHitInfo));
    fIChan = -1;
  }
 nextchan:
  if( NO_NEXT(fIChan, fNChan) )
    goto nextmod;
  Int_t chan = fEvData.GetNextChan(CRATE_SLOT(fHitInfo), fIChan);
  if( chan < fMod->lo or chan > fMod->hi )
    goto nextchan;  // Not one of my channels
  Int_t nhit = fEvData.GetNumHits(CRATE_SLOT(fHitInfo), chan);
  fHitInfo.chan = chan;
  fHitInfo.nhit = nhit;
  if( nhit == 0 ) {
    ostringstream ostr;
    ostr << "No hits on active "
         << (fHitInfo.type == Decoder::ChannelType::kADC ? "ADC" : "TDC")
         << " channel. Should never happen. Decoder bug. Call expert.";
    throw std::logic_error(msg(ostr.str().c_str()));
  }
  // If multiple hits on a TDC channel take the one earliest in time.
  // For a common-stop TDC, this is actually the last hit.
  fHitInfo.hit = (fHitInfo.type == Decoder::ChannelType::kCommonStopTDC) ? nhit - 1 : 0;
  // Determine default logical channel
  Int_t lchan = fMod->ConvertToLogicalChannel(chan);
  if( lchan < 0 or lchan >= size() ) {
    ostringstream ostr;
    ostr << "Illegal logical detector channel " << lchan << "."
         << "Must be between 0 and " << size() << ". Fix database.";
    throw std::invalid_argument(msg(ostr.str().c_str()));
  }
  fHitInfo.lchan = lchan;

  return *this;
}

//_____________________________________________________________________________
void THaDetMap::Iterator::reset()
{
  // Reset iterator to first element, if any

  fNMod = fDetMap.GetSize();
  fIMod = fIChan = -1;
  fHitInfo.reset();
  ++(*this);
}

//_____________________________________________________________________________
THaDetMap::Iterator& THaDetMap::MultiHitIterator::operator++()
{
  // Advance iterator to next active hit or channel, allowing multiple hits
  // per channel
 nexthit:
  if( LOOPDONE(fIHit, fHitInfo.nhit) ) {
    Iterator::operator++();
    if( !*this )
      return *this;
    fIHit = -1;
  }
  if( NO_NEXT(fIHit, fHitInfo.nhit) )
    goto nexthit;
  fHitInfo.hit = fIHit; // overwrite value of single-hit case (0 or nhit-1)

  return *this;
}

//_____________________________________________________________________________
ClassImp(THaDetMap)
