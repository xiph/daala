#!/usr/bin/env python
from collections import deque
import sys
import numpy as np
from skimage import color
import y4m

# Assuming BT.709
yuv2rgb = np.array([
    [1., 0., 1.28033], [1., -0.21482, -0.38059], [1., 2.12798, 0.]
])

# Simple box filter
box2 = np.ones((2, 2))


def decode_y4m_buffer(frame):
    W, H, buf = frame.headers['W'], frame.headers['H'], frame.buffer
    A, div2 = W * H, (H // 2, W // 2)
    Y = (np.ndarray((H, W), 'uint8', buf) - 16.) / 219.
    Cb = (np.ndarray(div2, 'uint8', buf, A) - 128.) / 224.
    Cr = (np.ndarray(div2, 'uint8', buf, A + A // 4) - 128.) / 224.
    YCbCr444 = np.dstack((Y, np.kron(Cb, box2), np.kron(Cr, box2)))
    return np.dot(YCbCr444, yuv2rgb.T)


scores = []


def process_pair(ref, recons):
    ref_lab = color.rgb2lab(decode_y4m_buffer(ref))
    recons_lab = color.rgb2lab(decode_y4m_buffer(recons))
    # "Color Image Quality Assessment Based on CIEDE2000"
    # Yang Yang, Jun Ming and Nenghai Yu, 2012
    # http://dx.doi.org/10.1155/2012/273723
    dE = color.deltaE_ciede2000(ref_lab, recons_lab, kL=0.65, kC=1.0, kH=4.0)
    scores.append(45. - 20. * np.log10(dE.mean()))
    print('%08d: %2.4f' % (ref.count, scores[-1]))


ref_frames = deque()
recons_frames = deque()


def process_ref(frame):
    ref_frames.append(frame)
    if recons_frames:
        process_pair(ref_frames.popleft(), recons_frames.popleft())


def process_recons(frame):
    recons_frames.append(frame)
    if ref_frames:
        process_pair(ref_frames.popleft(), recons_frames.popleft())


def main(args):
    OPENING = 'Opening %s...'
    BLOCK_SIZE = 4 * 1024 * 1024
    ref_parser = y4m.Reader(process_ref)
    recons_parser = y4m.Reader(process_recons)
    print(OPENING % args[1])
    with open(args[1], 'r') as ref:
        print(OPENING % args[2])
        with open(args[2], 'r') as recons:
            try:
                ref_buf, recons_buf = ref.buffer, recons.buffer
            except:
                ref_buf, recons_buf = ref, recons
            while True:
                data = ref_buf.read(BLOCK_SIZE)
                if not data: break
                ref_parser.decode(data)
                data = recons_buf.read(BLOCK_SIZE)
                if not data: break
                recons_parser.decode(data)
    print('Total: %2.4f' % np.array(scores).mean())


if __name__ == '__main__':
    main(sys.argv)
