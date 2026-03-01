/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    November 2013
  PURPOSE: A priority queue for PGM beliefs/distributions/messages.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "messagequeue.hpp"

// standard headers
#include <iostream>

using namespace std;

//-----------------------------------------------------------------------

MessageQueue::MsgPriority::MsgPriority(
  const Idx2& theFromTo,
  double theDelta)
    : fromTo(theFromTo), delta(theDelta) {}

void MessageQueue::MsgPriority::show() const {
  std::cout << fromTo << " " << delta << std::endl;
} // show

bool
MessageQueue::MsgCompare::operator() (
  const MsgPriority& msg1,
  const MsgPriority& msg2) const {
  return msg1.delta < msg2.delta;
} // operator()

unsigned
MessageQueue::size() {
  return msg2Delta.size();
} // size

void
MessageQueue::insert(
  const Idx2& idx,
  double val) {
  //cout << "Request to insert: " << idx << " " << val << endl;
  //cout << "Current queue size: " << size() << endl;
  auto mapItr = msg2Delta.find(idx);
  if ( mapItr != msg2Delta.end() ) { // existing message
    //cout << "Already exist with delta: " << mapItr->second->delta << "\n";
    val += mapItr->second->delta;
    //if ( val < mapItr->second->delta ) {return;} // do nothing
    //cout << "Erasing smaller existing\n";
    deltas.erase(mapItr->second); // erase existing smaller delta
    // dont need to erase the msg2Delta, will overwrite it
  } // if
  msg2Delta[idx] = deltas.insert( MsgPriority(idx,val) );
  //cout << "inserting " << idx << " " << val << ", new queue size: " << size() << endl;
  //show();
} //insert

double
MessageQueue::priority(
  const Idx2& idx){
  auto itr = msg2Delta.find(idx);
  if ( itr == msg2Delta.end() ) {return 0.0; } // if
  else {return itr->second->delta;} // else
} // priority

void
MessageQueue::remap(
  std::vector<unsigned> movedTo) {
  for (auto itr = msg2Delta.begin(); itr != msg2Delta.end(); itr++) {
    // get the info
    Idx2 oldIdx = itr->first;
    double val = itr->second->delta;
    // remove the old ones
    deltas.erase(itr->second);
    msg2Delta.erase(itr);
    // and add the new ones
    Idx2 newIdx(movedTo[oldIdx.first], movedTo[oldIdx.second]);
    msg2Delta[newIdx] = deltas.insert( MsgPriority(newIdx,val) );
  } // for
} // remap

MessageQueue::MsgPriority
MessageQueue::pop() {
  auto itr = --deltas.end();
  MsgPriority top = *itr;
  deltas.erase(itr);
  msg2Delta.erase( msg2Delta.find(top.fromTo) );
  return top;
} // pop

void
MessageQueue::clear() {
  msg2Delta.clear();
  deltas.clear();
} // clear

void
MessageQueue::show() const {
  for (auto msg : deltas) {
    msg.show();
  } // for
} // show
