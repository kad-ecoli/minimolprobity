// name: CliqueList.h
// author: J. Michael Word
// date written: 2/7/98
// purpose: Sequences containing collections of motion pointers
//          which refer to interacting residues, along with singletons

// **************************************************************
// NOTICE: This is free software and the source code is freely
// available. You are free to redistribute or modify under the
// conditions that (1) this notice is not removed or modified
// in any way and (2) any modified versions of the program are
// also available for free.
//               ** Absolutely no Warranty **
// Copyright (C) 1999 J. Michael Word
// **************************************************************

#ifndef CLIQUELIST_H
#define CLIQUELIST_H 1

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include "Mover.h"

class AtomPositions; // SJ 09/25/2015

class CliqueList {
public:
   CliqueList() {}
  ~CliqueList() {}

   CliqueList(const CliqueList& cl) : _cliques(cl._cliques),
                                      _singles(cl._singles) {}

   CliqueList& operator=(const CliqueList& cl) {
      if (this != &cl) {
	 _cliques = cl._cliques; _singles = cl._singles;
      }
      return *this;
   }

   void insertClique(const std::list<MoverPtr>& lst) { _cliques.push_front(lst); }
   void insertSingleton(const MoverPtr    ptr) { _singles.push_front(ptr); }

   void sortSingletonsByDescr();

   const std::list< std::list<MoverPtr> >& cliques() { return _cliques; }
   const std::list<MoverPtr>& singles() { return _singles; }

   int numCliques()    { return static_cast<int>(_cliques.size()); }
   int numSingletons() { return static_cast<int>(_singles.size()); }

   void describe(std::ostream& os) const;
   void formatSingles(std::vector<std::string>& cliqueNotes, AtomPositions& xyz) const; // SJ - 09/25/2015  added the last argument as that is need for doing the final flip
    void formatClique(std::vector<std::string>& cliqueNotes, int c, AtomPositions& xyz) const; // SJ - 09/25/2015  added the last argument as that is need for doing the final flip

private:
	std::list< std::list<MoverPtr> > _cliques;
	std::list<MoverPtr>        _singles;
};
#endif
