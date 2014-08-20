/**
 * @file neighbor_search_rules.hpp
 * @author Ryan Curtin
 *
 * Defines the pruning rules and base case rules necessary to perform a
 * tree-based search (with an arbitrary tree) for the NeighborSearch class.
 *
 * This file is part of MLPACK 1.0.9.
 *
 * MLPACK is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * MLPACK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details (LICENSE.txt).
 *
 * You should have received a copy of the GNU General Public License along with
 * MLPACK.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef __MLPACK_METHODS_NEIGHBOR_SEARCH_NEIGHBOR_SEARCH_RULES_HPP
#define __MLPACK_METHODS_NEIGHBOR_SEARCH_NEIGHBOR_SEARCH_RULES_HPP

#include "ns_traversal_info.hpp"

namespace mlpack {
namespace neighbor {

template<typename SortPolicy, typename MetricType, typename TreeType>
class NeighborSearchRules
{
 public:
  NeighborSearchRules(const typename TreeType::Mat& referenceSet,
                      const typename TreeType::Mat& querySet,
                      arma::Mat<long>& neighbors,
                      arma::mat& distances,
                      MetricType& metric);
  /**
   * Get the distance from the query point to the reference point.
   * This will update the "neighbor" matrix with the new point if appropriate
   * and will track the number of base cases (number of points evaluated).
   *
   * @param queryIndex Index of query point.
   * @param referenceIndex Index of reference point.
   */
  double BaseCase(const long queryIndex, const long referenceIndex);

  /**
   * Get the score for recursion order.  A low score indicates priority for
   * recursion, while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).
   *
   * @param queryIndex Index of query point.
   * @param referenceNode Candidate node to be recursed into.
   */
  double Score(const long queryIndex, TreeType& referenceNode);

  /**
   * Re-evaluate the score for recursion order.  A low score indicates priority
   * for recursion, while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).  This is used when the score has already
   * been calculated, but another recursion may have modified the bounds for
   * pruning.  So the old score is checked against the new pruning bound.
   *
   * @param queryIndex Index of query point.
   * @param referenceNode Candidate node to be recursed into.
   * @param oldScore Old score produced by Score() (or Rescore()).
   */
  double Rescore(const long queryIndex,
                 TreeType& referenceNode,
                 const double oldScore) const;

  /**
   * Get the score for recursion order.  A low score indicates priority for
   * recursionm while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).
   *
   * @param queryNode Candidate query node to recurse into.
   * @param referenceNode Candidate reference node to recurse into.
   */
  double Score(TreeType& queryNode, TreeType& referenceNode);

  /**
   * Re-evaluate the score for recursion order.  A low score indicates priority
   * for recursion, while DBL_MAX indicates that the node should not be recursed
   * into at all (it should be pruned).  This is used when the score has already
   * been calculated, but another recursion may have modified the bounds for
   * pruning.  So the old score is checked against the new pruning bound.
   *
   * @param queryNode Candidate query node to recurse into.
   * @param referenceNode Candidate reference node to recurse into.
   * @param oldScore Old score produced by Socre() (or Rescore()).
   */
  double Rescore(TreeType& queryNode,
                 TreeType& referenceNode,
                 const double oldScore) const;

  //! Get the number of base cases that have been performed.
  long BaseCases() const { return baseCases; }
  //! Modify the number of base cases that have been performed.
  long& BaseCases() { return baseCases; }

  //! Get the number of scores that have been performed.
  long Scores() const { return scores; }
  //! Modify the number of scores that have been performed.
  long& Scores() { return scores; }

  //! Convenience typedef.
  typedef NeighborSearchTraversalInfo<TreeType> TraversalInfoType;

  //! Get the traversal info.
  const TraversalInfoType& TraversalInfo() const { return traversalInfo; }
  //! Modify the traversal info.
  TraversalInfoType& TraversalInfo() { return traversalInfo; }

 private:
  //! The reference set.
  const typename TreeType::Mat& referenceSet;

  //! The query set.
  const typename TreeType::Mat& querySet;

  //! The matrix the resultant neighbor indices should be stored in.
  arma::Mat<long>& neighbors;

  //! The matrix the resultant neighbor distances should be stored in.
  arma::mat& distances;

  //! The instantiated metric.
  MetricType& metric;

  //! The last query point BaseCase() was called with.
  long lastQueryIndex;
  //! The last reference point BaseCase() was called with.
  long lastReferenceIndex;
  //! The last base case result.
  double lastBaseCase;

  //! The number of base cases that have been performed.
  long baseCases;
  //! The number of scores that have been performed.
  long scores;

  //! Traversal info for the parent combination; this is updated by the
  //! traversal before each call to Score().
  TraversalInfoType traversalInfo;

  /**
   * Recalculate the bound for a given query node.
   */
  double CalculateBound(TreeType& queryNode) const;

  /**
   * Insert a point into the neighbors and distances matrices; this is a helper
   * function.
   *
   * @param queryIndex Index of point whose neighbors we are inserting into.
   * @param pos Position in list to insert into.
   * @param neighbor Index of reference point which is being inserted.
   * @param distance Distance from query point to reference point.
   */
  void InsertNeighbor(const long queryIndex,
                      const long pos,
                      const long neighbor,
                      const double distance);
};

}; // namespace neighbor
}; // namespace mlpack

// Include implementation.
#include "neighbor_search_rules_impl.hpp"

#endif // __MLPACK_METHODS_NEIGHBOR_SEARCH_NEIGHBOR_SEARCH_RULES_HPP
