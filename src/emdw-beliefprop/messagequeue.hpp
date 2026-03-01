#ifndef MESSAGEQUEUE_HPP
#define MESSAGEQUEUE_HPP

/*******************************************************************************

  AUTHOR:  JA du Preez
  DATE:    November 2013
  PURPOSE: A priority queue for PGM beliefs/distributions/messages.

  (Copyright University of Stellenbosch, all rights reserved.)

*******************************************************************************/

// emdw headers
#include "emdw.hpp"

// standard headers
#include <map>
#include <set>

/******************************************************************************/

/**
 * @brief A priority queue for messages / distributions where we can
 * also change the priority/delta value of a message while it still is
 * in the queue.
 *
 * A message/distribution in the queue consists of an index pair (an
 * Idx2 i.e. a pair<unsigned,unsigned>) indicating from which factor
 * to which factor the message was sent, and a (positive) priority
 * value (typically indicating the degree of change relative to an
 * earlier version of this message).
 *
 * NOTE: If the need arises to prioritise on different criteria than
 * positive bigger is better, we might want to templatize on the
 * comparison functor.
 */
class MessageQueue{

public:

  class MsgPriority;
  class MsgCompare;

  /**
   * @brief the number of entries in the buffer - check this before
   * popping.
   */
  unsigned
  size();

  /**
   * @brief Inserts a message into the buffer - will update an
   * existing one if the new one has bigger delta.
   */
  void
  insert(
    const Idx2& idx,
    double val);

  /**
   * returns the priority value if in queue, 0 otherwise
   */
  double
  priority(
    const Idx2& idx);

  /**
   * @brief relocates a particular message to new indices
   *
   * @param movedTo The index that a specific factor has relocated to.
   */
  void
  remap(
    std::vector<unsigned> movedTo);

  /**
   * @brief pops and returns the top scoring message from the priority buffer
   */
  MsgPriority
  pop();

  /**
   * @brief erases the buffer.
   */
  void
  clear();

  // just for debugging
  void
  show() const;

  /// @brief helper class to encapsulate msg index and priority.
  class MsgPriority{
  public:

    MsgPriority(
      const Idx2& theFromTo,
      double theDelta);

    void
    show() const;

    Idx2 fromTo;
    double delta;
  }; // MsgPriority

  /// @brief helper class to compare messages on priority.
  class MsgCompare { // compares two messages to find smaller
  public:
    // returns true if msg1.delta < msg2.delta
    bool
    operator() (
      const MsgPriority& msg1,
      const MsgPriority& msg2) const;
  }; // MsgCompare

private:

  typedef std::multiset<MsgPriority,MsgCompare>::iterator Itr;

  // This points to the messages inside our priority queue. We use it
  // to access/update those messages, thereby re-arranging the queue.
  std::map<Idx2,Itr> msg2Delta;
  // This one serves as a priority queue, but one where we can reach
  // inside and update its elements with new information
  std::multiset<MsgPriority,MsgCompare> deltas;

}; // MessageQueue

#endif // MESSAGEQUEUE_HPP
