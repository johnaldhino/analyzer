#ifndef Podd_THaVDCPointPair_h_
#define Podd_THaVDCPointPair_h_

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// THaVDCPointPair                                                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TObject.h"
#include "THaVDCCluster.h"   // for chi2_t

class THaVDCPoint;
class THaTrack;

class THaVDCPointPair : public TObject {

public:
  THaVDCPointPair( THaVDCPoint* lp, THaVDCPoint* up, Double_t spacing, Double_t XRes, Double_t YRes)
    : fLowerPoint(lp), fUpperPoint(up), fSpacing(spacing), fXRes(XRes), fYRes(YRes), fError(1e38),
      fStatus(0) {}
  virtual ~THaVDCPointPair() {}

  void            Analyze();
  void            Associate( THaTrack* track );
  VDC::chi2_t     CalcChi2() const;
  virtual Int_t   Compare( const TObject* ) const;
  Double_t        GetError()   const { return fError; }
  THaVDCPoint*    GetLower()   const { return fLowerPoint; }
  THaVDCPoint*    GetUpper()   const { return fUpperPoint; }
  Double_t        GetSpacing() const { return fSpacing; }
  Int_t           GetStatus()  const { return fStatus; }
  THaTrack*       GetTrack()   const;
  Bool_t          HasUsedCluster() const;
  virtual Bool_t  IsSortable() const { return kTRUE; }
  virtual void    Print( Option_t* opt="" ) const;
  void            Release();
  void            SetStatus( Int_t i ) { fStatus = i; }
  void            Use();


  Double_t CalcError();

  
  static Double_t CalcErrorEst( THaVDCPoint* lowerPoint,
			     THaVDCPoint* upperPoint,
			     Double_t spacing );
  

  static void CalcXY( THaVDCPoint* lowerPoint,
			   THaVDCPoint* upperPoint,
			   Double_t spacing, Double_t &UV12X, Double_t &UV12Y, Double_t &UV12PX, Double_t &UV12PY, Double_t &UV21X, Double_t &UV21Y, Double_t &UV21PX, Double_t &UV21PY);


  static Double_t GetProjectedDistance( THaVDCPoint* here,
					THaVDCPoint* there,
					Double_t spacing, Double_t XRes = 0.009, Double_t YRes = 0.006); // default values for x and y resolution here established from data


  static void GetProjectedXY( THaVDCPoint* here,
			      THaVDCPoint* there,
			      Double_t spacing,
			      Double_t &X,
			      Double_t &Y,
			      Double_t &PX,
			      Double_t &PY);

  static Double_t GetTimeOffsetDifference(THaVDCPoint* here,
					  THaVDCPoint* there);

protected:

  THaVDCPoint*    fLowerPoint;  // Lower UV point
  THaVDCPoint*    fUpperPoint;  // Upper UV point
  Double_t        fSpacing;     // Spacing between lower and upper chambers [m]
  Double_t        fXRes;        // Resolution of X-PX [m]
  Double_t        fYRes;        // Resolution of Y-PY [m]
  Double_t        fError;       // Goodness of match between the points
  Int_t           fStatus;      // Status flag

private:
  THaVDCPointPair();

  ClassDef(THaVDCPointPair,0)     // Pair of lower/upper VDC points
};

//////////////////////////////////////////////////////////////////////////////

#endif
