// name: AtomConn.h
// author: J. Michael Word
// date written: 7/15/97
// purpose: Connected atoms and proton placement plans

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef ATOMCONN_H
#define ATOMCONN_H 1

#include <vector>
#include "Point3d.h"
#include "ElementInfo.h"
#include "utility.h"

// -----------------
//  atom plan flags
// -----------------

//NOTE:  STRICTALTFLAG is used only to prevent alpha H from worrying about sidechain alternate conformations
#define BONDBUMPFLAG    (1<<0)
#define STRICTALTFLAG   (1<<1)
#define   O2PRIMEFLAG   (1<<2)
#define NOO2PRIMEFLAG   (1<<3)
#define XTRAFLAG        (1<<4)
#define UNSUREDROPFLAG  (1<<5)
#define ROTATEFLAG      (1<<6)
#define ROTATEONDEMAND  (1<<7)
#define AROMATICFLAG    (1<<8)
#define HBACCEPTORFLAG  (1<<9)
#define HBDONORFLAG    (1<<10)
#define IFNOPO4        (1<<11)
#define NH3FLAG        (1<<12)
#define ISACOFLAG      (1<<13)

#define NEGCHARGEFLAG  (1<<14)
#define POSCHARGEFLAG  (1<<15)

// the following two flags are used to mark polar hydrogens where
// XPLOR has a different naming convention

#define NOTXPLORNAME    (1<<16)
#define    XPLORNAME    (1<<17)

// the following two flags are used to mark new vs. old naming conventions
#define USEOLDNAMES     (1<<18)
#define USENEWNAMES     (1<<19)

// the following two flags are used to allow two hydrogens 
// to be built on Calpha of backbone models
#define    NOTBBMODEL	(1<<20)
#define BACKBONEMODEL	(1<<21)

// -----------------------------------------
//  Keep track of the atoms connected to an atom 
// -----------------------------------------
class AtomConn {
public:
	AtomConn(const std::string& name, int ord) : _name(name), _order(ord) {}
  virtual ~AtomConn() {}

   AtomConn(const AtomConn& a)
	   : _name(a._name), _order(a._order) { _neighbor = a._neighbor; }
   AtomConn& operator=(const AtomConn& a);

   const std::string& name() const { return _name; }
   int num_conn() const { return static_cast<int>(_neighbor.size()); }
   int order() const { return _order; }

   void addConn(const std::string& c) { _neighbor.push_back(c); }

   void limitConnections(int n);
  
   // return the "i"th connection [0..num_conn-1] for this atom or empty string if out of range
   const std::string& conn(int i) const;
   
   // relational operator allow sequences of AtomConns to be sorted
   bool operator<(const AtomConn& a) const { return _order < a._order; }

private:

   std::string       _name;      // name of this atom
   int          _order;     // number used when sorting AtomConns and also used by algorithms as an index
   std::vector<std::string> _neighbor;  // list of connected atom names
};

// ---------------------------
//  new atom construction map
// ---------------------------
class atomPlacementPlan : public AtomConn {
public:
  /// @brief Construct an atom placement plan
  /// @param [in] Connection type (See StdRedH.h: 1=HXR3, 2=H2XR2, 3=H3XR, 4-5=HXR2, 6=HXY)
  /// @param [in] e ElementInfo description of the type of atom to be placed (Hydrogen).
  /// @param [in] c Atom connection structure to base the plan on
  /// @param [in] d Proton to heavy atom distance
  /// @param [in] a1 First angle to consider during placement
  /// @param [in] a2 Second angle to consider during placement
  /// @param [in] f Flags (see atom plan flags at the top of this file)
  atomPlacementPlan(int t, ElementInfo& e, const AtomConn &c,
			 float d, float a1, float a2, int f)
     : AtomConn(c), _type(t), _elem(e),
                    _dist(d), _ang1(a1), _ang2(a2), _flags(f) {}

  virtual ~atomPlacementPlan() {}

  /// @brief Copy constructor
  atomPlacementPlan(const atomPlacementPlan& a);
  atomPlacementPlan& operator=(const atomPlacementPlan& a);

  int type() const { return _type;  }
  const ElementInfo& elem() const { return _elem;  }
  float dist() const { return _dist;  }
  int hasFeature(int f) const { return _flags & f; }
  int addFeature(int f)       { _flags = _flags | f; return _flags; } // add feature for aromatic ring - Aram 07/18/12

  /// @brief Used by genHydrogens to attempt to place the hydrogen according to the plan
  /// @param [in] loc Locations of the heavy atoms to attempt placement with respect to.
  ///         The number of atoms required depends on the type of placement that is
  ///         being done.
  /// @param [out] hpos Location to place the hydrogen
  /// @return true on success, false on failure
  bool placeH(const std::vector<Point3d>& loc, Point3d& hpos) const;

  /// @brief Used by FlipMemo::finalize() to generate alternate positions for flipped H atoms.
  ///
  /// This method calls one of the private placement methods based on the locType,
  /// passing in a subset of the arguments and returning the results.  Not all of
  /// the parameters are used for each locType.
  ///   Note that this is a static method and does not require a constructed atomPlacementPlan.
  ///   There is a lot of magic going on here in terms of knowing which parameters to
  /// pass based on locType.
  /// @param [in] locType Connection type (See StdRedH.h: 1=HXR3, 2=H2XR2, 3=H3XR, 4-5=HXR2, 6=HXY)
  /// @param [in] a1 Position of the first atom to take into account during placement
  /// @param [in] a2 Position of the second atom to take into account during placement
  /// @param [in] a3 Position of the third atom to take into account during placement
  /// @param [in] a4 Position of the fourth atom to take into account during placement
  /// @param [in] len Distance from the heavy atom to the proton
  /// @param [in] ang Angle to take into account during placement
  /// @param [in] xtra Another parameter, whose meaning depends on locType
  /// @return Location of the hydrogen on success, (-999.9,-999.9,-999.9) on invalid type failure.
  static Point3d calcLoc(int locType,
                    const Point3d& a1, const Point3d& a2,
                    const Point3d& a3, const Point3d& a4,
                    float len, float ang, float xtra);

private:
  // 1: HXR3 - requires just 4 atom centers
  static Point3d type1position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    const Point3d& a3pos, const Point3d& a4pos,
                    float bondlen);
  // 2: H2XR2- three atoms and an angle
  static Point3d type2position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    const Point3d& a3pos,
                    float bondlen, float angle, float fudge);
  // 3: H3XR - three atoms an angle and dihedral
  static Point3d type3position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    const Point3d& a3pos,
                    float bondlen, float theta, float phi);
  // 4: HXR2 - three atoms and a fudge factor
  static Point3d type4position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    const Point3d& a3pos,
                    float bondlen, float fudge);
  // 5: HXR2 - three atoms and a fraction
  static Point3d type5position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    const Point3d& a3pos,
                    float bondlen, float fract);
  // 6: HXY  - (linear) just two atoms
  static Point3d type6position(
                    const Point3d& a1pos, const Point3d& a2pos,
                    float bondlen);

  int         _type;  // connection type
  ElementInfo _elem;  // type of atom
  float       _dist;  // proton to heavy atom distance
  float       _ang1;  // angle one
  float       _ang2;  // angle two
  int         _flags; // flags
};

#endif
