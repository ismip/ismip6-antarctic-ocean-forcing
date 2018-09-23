def get_res(config):
    res = config.getfloat('grid', 'dx')/1000.
    if res == int(res):
        res = int(res)
    res = '{}km'.format(res)
    return res
