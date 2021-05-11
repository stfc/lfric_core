##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################


'''PSyclone transformation script for the Dynamo0p3 API to apply
colouring, OpenMP and redundant computation to the level1 halo for
setval_* generically.

'''
from psyclone.transformations import Dynamo0p3ColourTrans, \
                                     Dynamo0p3OMPLoopTrans, \
                                     OMPParallelTrans, \
                                     Dynamo0p3RedundantComputationTrans

from psyclone.domain.lfric.function_space import FunctionSpace


def trans(psy):
    '''Applies PSyclone colouring, OpenMP and redundant computation
    transformations.

    '''
    ctrans = Dynamo0p3ColourTrans()
    otrans = Dynamo0p3OMPLoopTrans()
    oregtrans = OMPParallelTrans()
    rtrans = Dynamo0p3RedundantComputationTrans()

    setval_count = 0
    # Loop over all of the Invokes in the PSy object
    for invoke in psy.invokes.invoke_list:

        print("Transforming invoke '{0}' ...".format(invoke.name))
        schedule = invoke.schedule

        # Make setval_* compute redundantly to the level 1 halo if it
        # is in its own loop.
        for loop in schedule.loops():
            if loop.iteration_space == "dof":
                if len(loop.kernels()) != 1:
                    raise Exception(
                        "Expecting loop to contain 1 call but found '{0}'".
                        format(len(loop.kernels())))
                if loop.kernels()[0].name in ["setval_c", "setval_x"]:
                    setval_count += 1
                    schedule, _ = rtrans.apply(loop, options={"depth": 1})

        # Colour loops over cells unless they are on discontinuous
        # spaces or over dofs
        for loop in schedule.loops():
            if loop.iteration_space == "cell_column" \
                and loop.field_space.orig_name \
                    not in FunctionSpace.VALID_DISCONTINUOUS_NAMES:
                schedule, _ = ctrans.apply(loop)

        # Add OpenMP to loops unless they are over colours
        for loop in schedule.loops():
            if loop.loop_type != "colours":
                schedule, _ = oregtrans.apply(loop)
                schedule, _ = otrans.apply(loop, options={"reprod": True})

        # Take a look at what we've done
        print("Found {0} setval calls".format(setval_count))
        schedule.view()

    return psy
