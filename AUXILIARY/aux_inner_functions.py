def is_lowest_level(local_iter, string_as_iterable, base_string):
        if isinstance(local_iter, base_string):
            if string_as_iterable and len(local_iter) > 1:
                return False
            elif string_as_iterable and len(local_iter) <= 1:
                return True
            else:
                return False
        elif is_dict_callable(local_iter):
            return False
        elif is_list_callable(local_iter):
            return False
        else:
            return True
def is_dict_callable(local_iter):
    try:
        list(local_iter.keys())
        return True
    except:
        return False
def is_list_callable(local_iter):
    try:
        iter(local_iter)
        return True
    except:
        return False