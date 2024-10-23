#===========================================================================
#
# Section 392 Navigation and Mission Design
#
# Copyright (C) 2018, California Institute of Technology.
# U.S. Government Sponsorship under NASA Contract NAS7-03001 is acknowledged.
#
#===========================================================================

""" translatScaler class."""

from __future__ import print_function

from builtins import object

__version__ = '1.0.0'
__author__  = 'Juan Senent (392P) & Etienne Pellegrini (392M)'

#===========================================================================
# Place all imports after here.
#
from mmath import fabs

import Monte as M
#
# Place all imports before here.
#===========================================================================


#===========================================================================
class translatScaler(object):
   """ translatScaler class.

   This class defines a scaler for Cosmic using Juan Senent's strategy and
   inspired by the poincare correctorScaler. The constraint scales are computed using
   minimum bounds and user-defined tolerances.
   """

   #-----------------------------------------------------------------------
   def __init__(self, userTol={}, constraintScaleFactor=1.0E-5, defaultTol=1.0E-6):
      """ Constructor.

      = INPUT VARIABLES
      - userTol                dict: User-defined tolerances for user-defined constraints
      - constraintScaleFactor  float: The constraints scale factor
      """
      self.constraintScaleFactor = constraintScaleFactor
      self.userTol               = userTol
      self.defaultTol            = defaultTol

   #-----------------------------------------------------------------------
   def controlScale(self, control, refEpoch):
      """ Scales the control using the poincare corrector scaler strategy, which
      uses the maximum of the control bounds

      = INPUT VARIABLES
      - control   OptControl: The control to scale
      - refEpoch  Epoch: Not used
      """
      factor = max(fabs(control.minBound()), fabs(control.maxBound()))
      return 1.0 / factor

   #-----------------------------------------------------------------------
   def constraintScale(self, constraint):
      """ If the dictionary contains a user-defined tolerance for the constraint,
      use it. If not, use the absolute value of the minimum bound.

      = INPUT VARIABLES
      - constraint   OptConstraint: The constraint to scale
      """
      scales = []
      if constraint.name() in self.userTol.keys():
         scales.append(self.constraintScaleFactor / fabs(self.userTol[constraint.name()]))
      else:
         for (i, tol) in enumerate(constraint.minBounds()):
            if fabs(tol) > 0:
               scales.append(self.constraintScaleFactor / max(fabs(tol), fabs(constraint.maxBounds()[i])))
            else:
               # print('translatScaler: the constraint {0} does not have a '
                     # 'specified user tolerance or a minimum bound. Exiting.'.format(constraint.name()))
               # sys.exit()
               print('translatScaler: WARNING! the constraint {0} does not have a '
                     'specified user tolerance or a minimum bound. Using the default tol.'.format(constraint.name()))
               scales.append(M.UnitDbl(self.constraintScaleFactor / fabs(self.defaultTol), constraint.units()[0].reciprocal()))
               print('Scale: {0}'.format(scales[-1]))

      return scales

   #-----------------------------------------------------------------------
   def costScale(self, func):
      """ The cost is not scaled in this scaler. """
      return []

#===========================================================================
