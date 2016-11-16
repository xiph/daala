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
    W, H, C, buf = frame.headers['W'], frame.headers['H'], frame.headers['C'], frame.buffer
    A, div2, dtype, scale = W * H, (H // 2, W // 2), 'uint8', 1.
    if C.endswith('p10'):
        dtype, scale, A = 'uint16', 4., A * 2
    Y = (np.ndarray((H, W), dtype, buf) - 16. * scale) / (219. * scale)
    if C.startswith('420'):
        Cb = (np.ndarray(div2, dtype, buf, A) - 128. * scale) / (224. * scale)
        Cr = (np.ndarray(div2, dtype, buf, A + A // 4) - 128. * scale) / (224. * scale)
        YCbCr444 = np.dstack((Y, np.kron(Cb, box2), np.kron(Cr, box2)))
    else:
        Cb = (np.ndarray((H, W), dtype, buf, A) - 128. * scale) / (224. * scale)
        Cr = (np.ndarray((H, W), dtype, buf, A * 2) - 128. * scale) / (224. * scale)
        YCbCr444 = np.dstack((Y, Cb, Cr))
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

class Reader(y4m.Reader):
    def _frame_size(self):
        area = self._stream_headers['W'] * self._stream_headers['H']
        C = self._stream_headers['C']
        if C.startswith('420'):
            pixels = area * 3 // 2
        elif C.startswith('444'):
            pixels = area * 3
        else:
            raise Exception('Unknown pixel format: %s' % C)
        if self._stream_headers['C'].endswith('p10'):
            return 2 * pixels
        return pixels

def main(args):
    OPENING = 'Opening %s...'
    BLOCK_SIZE = 4 * 1024 * 1024
    ref_parser = Reader(process_ref)
    recons_parser = Reader(process_recons)
    print(OPENING % args[1])
    with open(args[1], 'r') as ref:
        print(OPENING % args[2])
        with open(args[2], 'r') as recons:
            try:
                ref_buf, recons_buf = ref.buffer, recons.buffer
            except:
                ref_buf, recons_buf = ref, recons
            while True:
                if not ref_frames:
                    data = ref_buf.read(BLOCK_SIZE)
                    if not data: break
                    ref_parser.decode(data)
                if not recons_frames:
                    data = recons_buf.read(BLOCK_SIZE)
                    if not data: break
                    recons_parser.decode(data)
    print('Total: %2.4f' % np.array(scores).mean())


if __name__ == '__main__':
    main(sys.argv)
