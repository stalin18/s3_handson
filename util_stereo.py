import math
import s3_handson.png as png
import struct

from s3_handson.util import *

class Writer:
    """
    PNG encoder in pure Python.
    """

    def __init__(self, width=None, height=None,
                 size=None,
                 greyscale=False,
                 alpha=False,
                 bitdepth=8,
                 palette=None,
                 transparent=None,
                 background=None,
                 gamma=None,
                 compression=None,
                 interlace=False,
                 bytes_per_sample=None, # deprecated
                 planes=None,
                 colormap=None,
                 maxval=None,
                 chunk_limit=2**20,
                 x_pixels_per_unit = None,
                 y_pixels_per_unit = None,
                 unit_is_meter = False):
        """
        Create a PNG encoder object.

        Arguments:

        width, height
          Image size in pixels, as two separate arguments.
        size
          Image size (w,h) in pixels, as single argument.
        greyscale
          Input data is greyscale, not RGB.
        alpha
          Input data has alpha channel (RGBA or LA).
        bitdepth
          Bit depth: from 1 to 16.
        palette
          Create a palette for a colour mapped image (colour type 3).
        transparent
          Specify a transparent colour (create a ``tRNS`` chunk).
        background
          Specify a default background colour (create a ``bKGD`` chunk).
        gamma
          Specify a gamma value (create a ``gAMA`` chunk).
        compression
          zlib compression level: 0 (none) to 9 (more compressed);
          default: -1 or None.
        interlace
          Create an interlaced image.
        chunk_limit
          Write multiple ``IDAT`` chunks to save memory.
        x_pixels_per_unit (pHYs chunk)
          Number of pixels a unit along the x axis
        y_pixels_per_unit (pHYs chunk)
          Number of pixels a unit along the y axis    
          With x_pixel_unit, give the pixel size ratio
        unit_is_meter (pHYs chunk)
          Indicates if unit is meter or not

        The image size (in pixels) can be specified either by using the
        `width` and `height` arguments, or with the single `size`
        argument.  If `size` is used it should be a pair (*width*,
        *height*).

        `greyscale` and `alpha` are booleans that specify whether
        an image is greyscale (or colour), and whether it has an
        alpha channel (or not).

        `bitdepth` specifies the bit depth of the source pixel values.
        Each source pixel value must be an integer between 0 and
        ``2**bitdepth-1``.  For example, 8-bit images have values
        between 0 and 255.  PNG only stores images with bit depths of
        1,2,4,8, or 16.  When `bitdepth` is not one of these values,
        the next highest valid bit depth is selected, and an ``sBIT``
        (significant bits) chunk is generated that specifies the
        original precision of the source image.  In this case the
        supplied pixel values will be rescaled to fit the range of
        the selected bit depth.

        The details of which bit depth / colour model combinations the
        PNG file format supports directly, are somewhat arcane
        (refer to the PNG specification for full details).  Briefly:
        "small" bit depths (1,2,4) are only allowed with greyscale and
        colour mapped images; colour mapped images cannot have bit depth
        16.

        For colour mapped images (in other words, when the `palette`
        argument is specified) the `bitdepth` argument must match one of
        the valid PNG bit depths: 1, 2, 4, or 8.  (It is valid to have a
        PNG image with a palette and an ``sBIT`` chunk, but the meaning
        is slightly different; it would be awkward to press the
        `bitdepth` argument into service for this.)

        The `palette` option, when specified, causes a colour mapped
        image to be created: the PNG colour type is set to 3; greyscale
        must not be set; alpha must not be set; transparent must not be
        set; the bit depth must be 1,2,4, or 8.  When a colour mapped
        image is created, the pixel values are palette indexes and
        the `bitdepth` argument specifies the size of these indexes
        (not the size of the colour values in the palette).

        The palette argument value should be a sequence of 3- or
        4-tuples.  3-tuples specify RGB palette entries; 4-tuples
        specify RGBA palette entries.  If both 4-tuples and 3-tuples
        appear in the sequence then all the 4-tuples must come
        before all the 3-tuples.  A ``PLTE`` chunk is created; if there
        are 4-tuples then a ``tRNS`` chunk is created as well.  The
        ``PLTE`` chunk will contain all the RGB triples in the same
        sequence; the ``tRNS`` chunk will contain the alpha channel for
        all the 4-tuples, in the same sequence.  Palette entries
        are always 8-bit.

        If specified, the `transparent` and `background` parameters must
        be a tuple with three integer values for red, green, blue, or
        a simple integer (or singleton tuple) for a greyscale image.

        If specified, the `gamma` parameter must be a positive number
        (generally, a float).  A ``gAMA`` chunk will be created.
        Note that this will not change the values of the pixels as
        they appear in the PNG file, they are assumed to have already
        been converted appropriately for the gamma specified.

        The `compression` argument specifies the compression level to
        be used by the ``zlib`` module.  Values from 1 to 9 specify
        compression, with 9 being "more compressed" (usually smaller
        and slower, but it doesn't always work out that way).  0 means
        no compression.  -1 and ``None`` both mean that the default
        level of compession will be picked by the ``zlib`` module
        (which is generally acceptable).

        If `interlace` is true then an interlaced image is created
        (using PNG's so far only interace method, *Adam7*).  This does
        not affect how the pixels should be presented to the encoder,
        rather it changes how they are arranged into the PNG file.
        On slow connexions interlaced images can be partially decoded
        by the browser to give a rough view of the image that is
        successively refined as more image data appears.

        .. note ::

          Enabling the `interlace` option requires the entire image
          to be processed in working memory.

        `chunk_limit` is used to limit the amount of memory used whilst
        compressing the image.  In order to avoid using large amounts of
        memory, multiple ``IDAT`` chunks may be created.
        """

        # At the moment the `planes` argument is ignored;
        # its purpose is to act as a dummy so that
        # ``Writer(x, y, **info)`` works, where `info` is a dictionary
        # returned by Reader.read and friends.
        # Ditto for `colormap`.

        width, height = check_sizes(size, width, height)
        del size

        if width <= 0 or height <= 0:
            raise ValueError("width and height must be greater than zero")
        if not isinteger(width) or not isinteger(height):
            raise ValueError("width and height must be integers")
        # http://www.w3.org/TR/PNG/#7Integers-and-byte-order
        if width > 2**32-1 or height > 2**32-1:
            raise ValueError("width and height cannot exceed 2**32-1")

        if alpha and transparent is not None:
            raise ValueError(
                "transparent colour not allowed with alpha channel")

        if bytes_per_sample is not None:
            warnings.warn('please use bitdepth instead of bytes_per_sample',
                          DeprecationWarning)
            if bytes_per_sample not in (0.125, 0.25, 0.5, 1, 2):
                raise ValueError(
                    "bytes per sample must be .125, .25, .5, 1, or 2")
            bitdepth = int(8*bytes_per_sample)
        del bytes_per_sample
        if not isinteger(bitdepth) or bitdepth < 1 or 16 < bitdepth:
            raise ValueError("bitdepth (%r) must be a positive integer <= 16" %
              bitdepth)

        self.rescale = None
        palette = check_palette(palette)
        if palette:
            if bitdepth not in (1,2,4,8):
                raise ValueError("with palette, bitdepth must be 1, 2, 4, or 8")
            if transparent is not None:
                raise ValueError("transparent and palette not compatible")
            if alpha:
                raise ValueError("alpha and palette not compatible")
            if greyscale:
                raise ValueError("greyscale and palette not compatible")
        else:
            # No palette, check for sBIT chunk generation.
            if alpha or not greyscale:
                if bitdepth not in (8,16):
                    targetbitdepth = (8,16)[bitdepth > 8]
                    self.rescale = (bitdepth, targetbitdepth)
                    bitdepth = targetbitdepth
                    del targetbitdepth
            else:
                assert greyscale
                assert not alpha
                if bitdepth not in (1,2,4,8,16):
                    if bitdepth > 8:
                        targetbitdepth = 16
                    elif bitdepth == 3:
                        targetbitdepth = 4
                    else:
                        assert bitdepth in (5,6,7)
                        targetbitdepth = 8
                    self.rescale = (bitdepth, targetbitdepth)
                    bitdepth = targetbitdepth
                    del targetbitdepth

        if bitdepth < 8 and (alpha or not greyscale and not palette):
            raise ValueError(
              "bitdepth < 8 only permitted with greyscale or palette")
        if bitdepth > 8 and palette:
            raise ValueError(
                "bit depth must be 8 or less for images with palette")

        transparent = check_color(transparent, greyscale, 'transparent')
        background = check_color(background, greyscale, 'background')

        # It's important that the true boolean values (greyscale, alpha,
        # colormap, interlace) are converted to bool because Iverson's
        # convention is relied upon later on.
        self.width = width
        self.height = height
        self.transparent = transparent
        self.background = background
        self.gamma = gamma
        self.greyscale = bool(greyscale)
        self.alpha = bool(alpha)
        self.colormap = bool(palette)
        self.bitdepth = int(bitdepth)
        self.compression = compression
        self.chunk_limit = chunk_limit
        self.interlace = bool(interlace)
        self.palette = palette
        self.x_pixels_per_unit = x_pixels_per_unit
        self.y_pixels_per_unit = y_pixels_per_unit
        self.unit_is_meter = bool(unit_is_meter)

        self.color_type = 4*self.alpha + 2*(not greyscale) + 1*self.colormap
        assert self.color_type in (0,2,3,4,6)

        self.color_planes = (3,1)[self.greyscale or self.colormap]
        self.planes = self.color_planes + self.alpha
        # :todo: fix for bitdepth < 8
        self.psize = (self.bitdepth/8) * self.planes

    def make_palette(self):
        """Create the byte sequences for a ``PLTE`` and if necessary a
        ``tRNS`` chunk.  Returned as a pair (*p*, *t*).  *t* will be
        ``None`` if no ``tRNS`` chunk is necessary.
        """

        p = array('B')
        t = array('B')

        for x in self.palette:
            p.extend(x[0:3])
            if len(x) > 3:
                t.append(x[3])
        p = tostring(p)
        t = tostring(t)
        if t:
            return p,t
        return p,None

    def write(self, outfile, rows):
        """Write a PNG image to the output file.  `rows` should be
        an iterable that yields each row in boxed row flat pixel
        format.  The rows should be the rows of the original image,
        so there should be ``self.height`` rows of ``self.width *
        self.planes`` values.  If `interlace` is specified (when
        creating the instance), then an interlaced PNG file will
        be written.  Supply the rows in the normal image order;
        the interlacing is carried out internally.

        .. note ::

          Interlacing will require the entire image to be in working
          memory.
        """

        if self.interlace:
            fmt = 'BH'[self.bitdepth > 8]
            a = array(fmt, itertools.chain(*rows))
            return self.write_array(outfile, a)

        nrows = self.write_passes(outfile, rows)
        if nrows != self.height:
            raise ValueError(
              "rows supplied (%d) does not match height (%d)" %
              (nrows, self.height))

    def write_passes(self, outfile, rows, packed=False):
        """
        Write a PNG image to the output file.

        Most users are expected to find the :meth:`write` or
        :meth:`write_array` method more convenient.
        
        The rows should be given to this method in the order that
        they appear in the output file.  For straightlaced images,
        this is the usual top to bottom ordering, but for interlaced
        images the rows should have already been interlaced before
        passing them to this function.

        `rows` should be an iterable that yields each row.  When
        `packed` is ``False`` the rows should be in boxed row flat pixel
        format; when `packed` is ``True`` each row should be a packed
        sequence of bytes.
        """

        # http://www.w3.org/TR/PNG/#5PNG-file-signature
        outfile.write(_signature)

        # http://www.w3.org/TR/PNG/#11IHDR
        write_chunk(outfile, 'IHDR',
                    struct.pack("!2I5B", self.width, self.height,
                                self.bitdepth, self.color_type,
                                0, 0, self.interlace))

        # See :chunk:order
        # http://www.w3.org/TR/PNG/#11gAMA
        if self.gamma is not None:
            write_chunk(outfile, 'gAMA',
                        struct.pack("!L", int(round(self.gamma*1e5))))

        # See :chunk:order
        # http://www.w3.org/TR/PNG/#11sBIT
        if self.rescale:
            write_chunk(outfile, 'sBIT',
                struct.pack('%dB' % self.planes,
                            *[self.rescale[0]]*self.planes))
        
        # :chunk:order: Without a palette (PLTE chunk), ordering is
        # relatively relaxed.  With one, gAMA chunk must precede PLTE
        # chunk which must precede tRNS and bKGD.
        # See http://www.w3.org/TR/PNG/#5ChunkOrdering
        if self.palette:
            p,t = self.make_palette()
            write_chunk(outfile, 'PLTE', p)
            if t:
                # tRNS chunk is optional. Only needed if palette entries
                # have alpha.
                write_chunk(outfile, 'tRNS', t)

        # http://www.w3.org/TR/PNG/#11tRNS
        if self.transparent is not None:
            if self.greyscale:
                write_chunk(outfile, 'tRNS',
                            struct.pack("!1H", *self.transparent))
            else:
                write_chunk(outfile, 'tRNS',
                            struct.pack("!3H", *self.transparent))

        # http://www.w3.org/TR/PNG/#11bKGD
        if self.background is not None:
            if self.greyscale:
                write_chunk(outfile, 'bKGD',
                            struct.pack("!1H", *self.background))
            else:
                write_chunk(outfile, 'bKGD',
                            struct.pack("!3H", *self.background))

        # http://www.w3.org/TR/PNG/#11pHYs
        if self.x_pixels_per_unit is not None and self.y_pixels_per_unit is not None:
            tup = (self.x_pixels_per_unit, self.y_pixels_per_unit, int(self.unit_is_meter))
            write_chunk(outfile, 'pHYs', struct.pack("!LLB",*tup))

        # http://www.w3.org/TR/PNG/#11IDAT
        if self.compression is not None:
            compressor = zlib.compressobj(self.compression)
        else:
            compressor = zlib.compressobj()

        # Choose an extend function based on the bitdepth.  The extend
        # function packs/decomposes the pixel values into bytes and
        # stuffs them onto the data array.
        data = array('B')
        if self.bitdepth == 8 or packed:
            extend = data.extend
        elif self.bitdepth == 16:
            # Decompose into bytes
            def extend(sl):
                fmt = '!%dH' % len(sl)
                data.extend(array('B', struct.pack(fmt, *sl)))
        else:
            # Pack into bytes
            assert self.bitdepth < 8
            # samples per byte
            spb = int(8/self.bitdepth)
            def extend(sl):
                a = array('B', sl)
                # Adding padding bytes so we can group into a whole
                # number of spb-tuples.
                l = float(len(a))
                extra = math.ceil(l / float(spb))*spb - l
                a.extend([0]*int(extra))
                # Pack into bytes
                l = group(a, spb)
                l = map(lambda e: reduce(lambda x,y:
                                           (x << self.bitdepth) + y, e), l)
                data.extend(l)
        if self.rescale:
            oldextend = extend
            factor = \
              float(2**self.rescale[1]-1) / float(2**self.rescale[0]-1)
            def extend(sl):
                oldextend(map(lambda x: int(round(factor*x)), sl))

        # Build the first row, testing mostly to see if we need to
        # changed the extend function to cope with NumPy integer types
        # (they cause our ordinary definition of extend to fail, so we
        # wrap it).  See
        # http://code.google.com/p/pypng/issues/detail?id=44
        enumrows = enumerate(rows)
        del rows

        # First row's filter type.
        data.append(0)
        # :todo: Certain exceptions in the call to ``.next()`` or the
        # following try would indicate no row data supplied.
        # Should catch.
        if sys.version_info[0] == 2:
            i,row = enumrows.next()
        else:
            i,row = next(enumrows)
        try:
            # If this fails...
            extend(row)
        except:
            # ... try a version that converts the values to int first.
            # Not only does this work for the (slightly broken) NumPy
            # types, there are probably lots of other, unknown, "nearly"
            # int types it works for.
            def wrapmapint(f):
                return lambda sl: f(map(int, sl))
            extend = wrapmapint(extend)
            del wrapmapint
            extend(row)

        for i,row in enumrows:
            # Add "None" filter type.  Currently, it's essential that
            # this filter type be used for every scanline as we do not
            # mark the first row of a reduced pass image; that means we
            # could accidentally compute the wrong filtered scanline if
            # we used "up", "average", or "paeth" on such a line.
            data.append(0)
            extend(row)
            if len(data) > self.chunk_limit:
                compressed = compressor.compress(tostring(data))
                if len(compressed):
                    write_chunk(outfile, 'IDAT', compressed)
                # Because of our very witty definition of ``extend``,
                # above, we must re-use the same ``data`` object.  Hence
                # we use ``del`` to empty this one, rather than create a
                # fresh one (which would be my natural FP instinct).
                del data[:]
        if len(data):
            compressed = compressor.compress(tostring(data))
        else:
            compressed = strtobytes('')
        flushed = compressor.flush()
        if len(compressed) or len(flushed):
            write_chunk(outfile, 'IDAT', compressed + flushed)
        # http://www.w3.org/TR/PNG/#11IEND
        write_chunk(outfile, 'IEND')
        return i+1

    def write_array(self, outfile, pixels):
        """
        Write an array in flat row flat pixel format as a PNG file on
        the output file.  See also :meth:`write` method.
        """

        if self.interlace:
            self.write_passes(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write_passes(outfile, self.array_scanlines(pixels))

    def write_packed(self, outfile, rows):
        """
        Write PNG file to `outfile`.  The pixel data comes from `rows`
        which should be in boxed row packed format.  Each row should be
        a sequence of packed bytes.

        Technically, this method does work for interlaced images but it
        is best avoided.  For interlaced images, the rows should be
        presented in the order that they appear in the file.

        This method should not be used when the source image bit depth
        is not one naturally supported by PNG; the bit depth should be
        1, 2, 4, 8, or 16.
        """

        if self.rescale:
            raise Error("write_packed method not suitable for bit depth %d" %
              self.rescale[0])
        return self.write_passes(outfile, rows, packed=True)

    def convert_pnm(self, infile, outfile):
        """
        Convert a PNM file containing raw pixel data into a PNG file
        with the parameters set in the writer object.  Works for
        (binary) PGM, PPM, and PAM formats.
        """

        if self.interlace:
            pixels = array('B')
            pixels.fromfile(infile,
                            (self.bitdepth/8) * self.color_planes *
                            self.width * self.height)
            self.write_passes(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write_passes(outfile, self.file_scanlines(infile))

    def convert_ppm_and_pgm(self, ppmfile, pgmfile, outfile):
        """
        Convert a PPM and PGM file containing raw pixel data into a
        PNG outfile with the parameters set in the writer object.
        """
        pixels = array('B')
        pixels.fromfile(ppmfile,
                        (self.bitdepth/8) * self.color_planes *
                        self.width * self.height)
        apixels = array('B')
        apixels.fromfile(pgmfile,
                         (self.bitdepth/8) *
                         self.width * self.height)
        pixels = interleave_planes(pixels, apixels,
                                   (self.bitdepth/8) * self.color_planes,
                                   (self.bitdepth/8))
        if self.interlace:
            self.write_passes(outfile, self.array_scanlines_interlace(pixels))
        else:
            self.write_passes(outfile, self.array_scanlines(pixels))

    def file_scanlines(self, infile):
        """
        Generates boxed rows in flat pixel format, from the input file
        `infile`.  It assumes that the input file is in a "Netpbm-like"
        binary format, and is positioned at the beginning of the first
        pixel.  The number of pixels to read is taken from the image
        dimensions (`width`, `height`, `planes`) and the number of bytes
        per value is implied by the image `bitdepth`.
        """

        # Values per row
        vpr = self.width * self.planes
        row_bytes = vpr
        if self.bitdepth > 8:
            assert self.bitdepth == 16
            row_bytes *= 2
            fmt = '>%dH' % vpr
            def line():
                return array('H', struct.unpack(fmt, infile.read(row_bytes)))
        else:
            def line():
                scanline = array('B', infile.read(row_bytes))
                return scanline
        for y in range(self.height):
            yield line()

    def array_scanlines(self, pixels):
        """
        Generates boxed rows (flat pixels) from flat rows (flat pixels)
        in an array.
        """

        # Values per row
        vpr = self.width * self.planes
        stop = 0
        for y in range(self.height):
            start = stop
            stop = start + vpr
            yield pixels[start:stop]

    def array_scanlines_interlace(self, pixels):
        """
        Generator for interlaced scanlines from an array.  `pixels` is
        the full source image in flat row flat pixel format.  The
        generator yields each scanline of the reduced passes in turn, in
        boxed row flat pixel format.
        """

        # http://www.w3.org/TR/PNG/#8InterlaceMethods
        # Array type.
        fmt = 'BH'[self.bitdepth > 8]
        # Value per row
        vpr = self.width * self.planes
        for xstart, ystart, xstep, ystep in _adam7:
            if xstart >= self.width:
                continue
            # Pixels per row (of reduced image)
            ppr = int(math.ceil((self.width-xstart)/float(xstep)))
            # number of values in reduced image row.
            row_len = ppr*self.planes
            for y in range(ystart, self.height, ystep):
                if xstep == 1:
                    offset = y * vpr
                    yield pixels[offset:offset+vpr]
                else:
                    row = array(fmt)
                    # There's no easier way to set the length of an array
                    row.extend(pixels[0:row_len])
                    offset = y * vpr + xstart * self.planes
                    end_offset = (y+1) * vpr
                    skip = self.planes * xstep
                    for i in range(self.planes):
                        row[i::self.planes] = \
                            pixels[offset+i:end_offset:skip]
                    yield row
					
# Reads a .txt file containing the reported runtime in seconds for a Middlebury
# result as a float value as text. Returns the time as a float.
def ReadMiddlebury2014TimeFile(path):
    time = -1
    with open(path, 'rb') as time_file:
        text = time_file.read().decode('UTF-8').strip()
        try:
            time = float(text)
        except ValueError:
            raise Exception('Cannot parse time file: ' + path)
    return time


# Returns a dict which maps the parameters to their values. The values (right
# side of the equal sign) are all returned as strings (and not parsed).
def ReadMiddlebury2014CalibFile(path):
    result = dict()
    with open(path, 'rb') as calib_file:
        for line in calib_file.readlines():
            line = line.decode('UTF-8').rstrip('\n')
            if len(line) == 0:
                continue
            eq_pos = line.find('=')
            if eq_pos < 0:
                raise Exception('Cannot parse Middlebury 2014 calib file: ' + path)
            result[line[:eq_pos]] = line[eq_pos + 1:]
    return result


# Writes a calib.txt file according to the Middlebury format, given the required
# values.
def WriteMiddlebury2014CalibFile(path,
                                 left_fx, left_fy, left_cx, left_cy,
                                 right_fx, right_fy, right_cx, right_cy,
                                 baseline_in_mm,
                                 width,
                                 height,
                                 ndisp):
    with open(path, 'wb') as calib_file:
        calib_file.write(StrToBytes('cam0=[' + str(left_fx) + ' 0 ' + str(left_cx) + '; 0 ' + str(left_fy) + ' ' + str(left_cy) + '; 0 0 1]\n'))
        calib_file.write(StrToBytes('cam1=[' + str(right_fx) + ' 0 ' + str(right_cx) + '; 0 ' + str(right_fy) + ' ' + str(right_cy) + '; 0 0 1]\n'))
        calib_file.write(StrToBytes('doffs=' + str(right_cx - left_cx) + '\n'))
        calib_file.write(StrToBytes('baseline=' + str(baseline_in_mm) + '\n'))
        calib_file.write(StrToBytes('width=' + str(width) + '\n'))
        calib_file.write(StrToBytes('height=' + str(height) + '\n'))
        calib_file.write(StrToBytes('ndisp=' + str(ndisp) + '\n'))


# Reads a .pfm file containing a disparity image in Middlebury format.
# Returns a 3-tuple (width, height, pixels), where pixels is a tuple of floats,
# ordered as in the PFM file (i.e., the bottommost row comes first).
def ReadMiddlebury2014PfmFile(path):
    with open(path, 'rb') as pfm_file:
        state = 0
        word = ''
        width = -1
        height = -1
        little_endian = True
        while True:
            character = pfm_file.read(1).decode('UTF-8')
            if not character:
                raise Exception('Cannot parse pfm file: unexpected end of file')
            elif character == '#' or character == ' ' or character == '\n' or character == '\r' or character == '\t':
                # Parse word
                if word != '':
                    if state == 0:
                        if word != 'Pf':
                            raise Exception('Cannot parse pfm file: header is not "Pf"')
                        state = 1
                    elif state == 1:
                        width = int(word)
                        state = 2
                    elif state == 2:
                        height = int(word)
                        state = 3
                    elif state == 3:
                        little_endian = float(word) < 0
                        break
                
                word = ''
                
                if character == '#':
                    # Skip comment.
                    pfm_file.readline()
                else:
                    # Skip whitespace
                    continue
            
            word += character
        
        # Read float buffer
        pixel_count = width * height
        endian_character = '<' if little_endian else '>'
        pixels = struct.unpack(endian_character + str(pixel_count) + 'f', pfm_file.read(4 * pixel_count))
    
    return (width, height, pixels)


# Writes a .pfm file containing a disparity image according to Middlebury format.
# Expects pixels as a list of floats
def WriteMiddlebury2014PfmFile(path, width, height, pixels):
    with open(path, 'wb') as pfm_file:
        pfm_file.write(StrToBytes('Pf\n'))
        pfm_file.write(StrToBytes(str(width) + ' ' + str(height) + '\n'))
        pfm_file.write(StrToBytes('-1\n'))  # negative number means little endian
        pfm_file.write(struct.pack('<' + str(len(pixels)) + 'f', *pixels))  # < means using little endian


# Converts a Middlebury .pfm disparity image to a Kitti .png disparity image.
def ConvertMiddlebury2014PfmToKitti2015Png(src_path, dest_path):
    (pfm_width, pfm_height, pfm_pixels) = ReadMiddlebury2014PfmFile(src_path)
    
    png_disp = []  # list of rows
    for y in range(pfm_height - 1, -1, -1):  # iterate in reverse order according to pfm format
        in_row = pfm_pixels[y * pfm_width : (y + 1) * pfm_width]
        out_row = []
        for value in in_row:
            if math.isinf(value):
                out_row.append(0)  # invalid value
            else:
                converted_value = max(1, int(round(256.0 * value)))
                if converted_value > 65535:
                    print('Warning: A disparity value of 256 or larger needed to be clamped in the conversion from PFM to Kitti PNG. File: ' + src_path)
                    converted_value = 255
                out_row.append(converted_value)
        png_disp.append(out_row)
    
    with open(dest_path, 'wb') as dest_png_file:
        png_writer = Writer(width=pfm_width, height=pfm_height, bitdepth=16, compression=9, greyscale=True)
        #png_writer = png.Writer(width=pfm_width, height=pfm_height, bitdepth=8, compression=9, greyscale=True)
        png_writer.write(dest_png_file, png_disp)


# Converts a Kitti .png disparity image to a Middlebury .pfm disparity image.
def ConvertKitti2015PngToMiddlebury2014Pfm(src_path, dest_path):
    # Read .png file.
    disp_reader = png.Reader(src_path)
    disp_data = disp_reader.read()
    if disp_data[3]['bitdepth'] != 16:
        raise Exception('bitdepth of ' + src_path + ' is not 16')
    
    width = disp_data[0]
    height = disp_data[1]
    
    # Get list of rows.
    disp_rows = list(disp_data[2])
    
    # Convert to Middlebury's PFM format.
    disp_float = []
    for y in range(len(disp_rows) - 1, -1, -1):  # iterate in reverse order according to pfm format
        input_line = disp_rows[y]
        for value in input_line:
            if value > 0:
                disp_float.append(float(value) / 256.0)
            else:
                disp_float.append(float('inf'))  # invalid value
    
    WriteMiddlebury2014PfmFile(dest_path, width, height, disp_float)
    
    disp_reader.close()
