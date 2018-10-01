
def get_horiz_res(config):
    hres = config.getfloat('grid', 'dx')/1000.
    if hres == int(hres):
        hres = int(hres)
    res = '{}km'.format(hres)
    return res


def get_res(config, extrap=True):
    hres = get_horiz_res(config)
    if extrap:
        vres = -config.getfloat('grid', 'dzExtrap')
    else:
        vres = -config.getfloat('grid', 'dzFinal')
    if vres == int(vres):
        vres = int(vres)
    res = '{}_x_{}m'.format(hres, vres)
    return res
