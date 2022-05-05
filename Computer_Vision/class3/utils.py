import numpy as np

def calc_r(cmax, g, b, diff):
    h = (60 * ((g - b) / diff) + 360) % 360
    s = (diff / cmax) * 100
    v = cmax * 100
    return h, s, v


def calc_g(r, cmax, b, diff):
    h = (60 * ((b - r) / diff) + 120) % 360
    s = (diff / cmax) * 100
    v = cmax * 100
    return h, s, v


def calc_b(r, g, cmax, diff):
    h = (60 * ((r - g) / diff) + 240) % 360
    s = (diff / cmax) * 100
    v = cmax * 100
    return h, s, v

def calc_h0(x, c):
    r0, g0, b0 = c, x, 0
    return r0, g0, b0


def calc_h1(x, c):
    r0, g0, b0 = x, c, 0
    return r0, g0, b0


def calc_h2(x, c):
    r0, g0, b0 = 0, c, x
    return r0, g0, b0


def calc_h3(x, c):
    r0, g0, b0 = 0, x, c
    return r0, g0, b0


def calc_h4(x, c):
    r0, g0, b0 = x, 0, c
    return r0, g0, b0


def calc_h5(x, c):
    r0, g0, b0 = c, 0, x
    return r0, g0, b0


rgb_to_calc = {0: calc_r, 1: calc_g, 2: calc_b}


h6_to_calc = {0: calc_h0, 1: calc_h1, 2: calc_h2, 3: calc_h3, 4: calc_h4, 5: calc_h5}

def rgb_to_hsv(img_arr):
    hsv_img = img_arr / 255
    cmaxes = np.amax(hsv_img, 2)
    cminis = np.amin(hsv_img, 2)
    diffs = cmaxes - cminis
    
    h, w, _ = hsv_img.shape
    for i in range(h):
        for j in range(w):
            if cmaxes[i, j] < 1e-4:
                h, s, v = 0, 0, 0
            elif diffs[i, j] < 1e-4:
                h, s, v = 0, 0, cmaxes[i, j] * 100
            else:
                max_idx = np.where(abs(hsv_img[i, j] - cmaxes[i, j]) < 1e-5)[0][0]
                h, s, v = rgb_to_calc[max_idx](*hsv_img[i, j], diffs[i, j])
            hsv_img[i, j] = np.array([h, s, v])
    return hsv_img
    
def hsv_to_rgb(img_arr):
    rgb_img = np.empty_like(img_arr)
    c_arr = img_arr[:, :, 1] * img_arr[:, :, 2] /1e4
    x_arr = c_arr * (1 - np.abs((img_arr[:, :, 0] / 60 % 2 - 1)))
    m_arr = (img_arr[:, :, 2] / 100) - c_arr

    h, w, _ = rgb_img.shape
    for i in range(h):
        for j in range(w):
            h = img_arr[i, j, 0]
            r0, g0, b0 = h6_to_calc[h//60](x_arr[i,j], c_arr[i,j])
            rgb_img[i, j] = np.array([r0, g0, b0])
    rgb_img = (rgb_img+np.dstack([m_arr]*3))*255
    np.clip(rgb_img, 0, 255)
    return rgb_img

if __name__ == '__main__':
    main()