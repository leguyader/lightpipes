# -*- coding: utf-8 -*-
"""
Separated backward_compat into it's own module to avoid circular
import dependencies (core needed misc.backward_compat, but misc needed core)

Created on Sun Jul 12 13:35:54 2020

@author: Leonard.Doyle
"""

import functools

from .field import Field

def backward_compatible(fn):
    """Decorator to wrap new-convention (field first) LP methods to also
    accept old convention (field last).
    This decorator should work for any new-convention functions that
    only changed the order of the field. If the order or number of other
    arguments has also changed, this decorator might fail and one has to
    manually take care of the order inside that function!
    
    If the user supplies only positional arguments (without keyword), the
    field is searched as first argument for new-convention functions. If not
    found, the field is searched as last argument (old convention) and swapped
    to front if found. If neither first nor last, an error is raised.
    
    The handling of all possible combinations of positional and keyword-
    supplied arguments is tricky, so it might fail in rare occasions/
    strange use cases (e.g. someone supplying kwargs even though all args
    are given and in the correct order).
    Therefore, if at least one argument is supplied via keyword and the first
    positional argument is not a field, an error is raised (swapping would be
    guessing). The solution is to either provide all positional, all keyword
    or update to new syntax putting Field as first argument.
    
    Example new-convention definition:
    @backward_compatible
    def CircScreen(Fin, R, x_shift=0, y_shift=0):
        #old convention was
        # CircScreen(R, x_shift, y_shift, Fin)
    
    Valid function calls then are:
        F = CircScreen(F, 1, 0, 0) #new style, all positional
        F = CircScreen(F, 1) #new style, all positional, some default value
        F = CircScreen(F, 1, y_shift=1) #new style, some default, some kwarg
        
        F = CircScreen(1, 0, 0, F) #old style, all positional
        F = CircScreen(R=1, x_shift=0, y_shift=0, Fin=F) #old, all kwargs
        
        F = CircScreen(R=1, x_shift=0, Fin=F, y_shift=0) #...
            # all kwargs, order irrelevant
    
    Invalid function calls (resolving cannot be guaranteed by decorator):
        F = CircScreen(1, 0, 0, Fin=F) #old style, at least one kwarg
        F = CircScreen(1, F, 0, 0) #plain wrong, neither old nor new style
    
    Grey area:
        F = CircScreen(1, F) #old style part positional, part default value
            #interestingly, this will not raise an error in decorator
            # but might fail in the function call since order not clear
    """
    @functools.wraps(fn)
    def fn_wrapper(*args, **kwargs):
        args = list(args) #make mutable
        if len(kwargs)==0:
            #no keyword args supplied, all args stricly by order/positional
            # -> easy, this decorator won't fail
            if not isinstance(args[0], Field):
                #first arg is not a field, either backward compat syntax or
                # complete usage error -> find out if Field is last, else error
                if isinstance(args[-1], Field):
                    #found field in last arg, push last arg to first:
                    args.insert(0, args.pop()) #[1,2,3,4] -> [4, 1, 2, 3]
                    #-> now all the variables contain what is expected in new style
                else:
                    raise ValueError(fn.__name__ + '(backward compatibility'
                                     + ' check): Field is neither first'
                                     + ' nor last parameter, please check '
                                     + 'syntax/usage.')
        else:
            if len(args)>0:
                #at least one argument is supplied by explicitly naming it,
                # while others are positiobal args.
                # In this case flipping the order or the rest seems
                # dangerous/wrong/unpredictable
                if isinstance(args[0], Field):
                    #all good, Field still first even though kwargs later
                    pass
                else:
                    raise ValueError(fn.__name__ + ' (backward comaptibility'
                                     + ' check): Field not first and kwargs '
                                     + 'used, please update to new '
                                     + 'syntax/usage.')
        
        return fn(*args, **kwargs)
    
    return fn_wrapper

