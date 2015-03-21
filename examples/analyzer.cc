#include <wx/wx.h>
#include <wx/dcbuffer.h>

extern "C" {
#include <daala/codec.h>
#include <daala/daaladec.h>
}

/*Smallest blocks are 4x4*/
# define OD_LOG_BSIZE0 (2)
/*There are 4 block sizes total (4x4, 8x8, 16x16, 32x32).*/
# define OD_NBSIZES    (4)
/*The maximum length of the side of a block.*/
# define OD_BSIZE_MAX  (1 << OD_LOG_BSIZE0 + OD_NBSIZES - 1)

# define OD_MAXI(a, b) ((a) ^ (((a) ^ (b)) & -((b) > (a))))
# define OD_MINI(a, b) ((a) ^ (((b) ^ (a)) & -((b) < (a))))
# define OD_CLAMPI(a, b, c) (OD_MAXI(a, OD_MINI(b, c)))

# define OD_SIGNMASK(a) (-((a) < 0))
# define OD_FLIPSIGNI(a, b) (((a) + OD_SIGNMASK(b)) ^ OD_SIGNMASK(b))
# define OD_DIV_ROUND(x, y) (((x) + OD_FLIPSIGNI((y) >> 1, x))/(y))

# define OD_BLOCK_SIZE4x4(bsize, bstride, bx, by) \
 ((bsize)[((by) >> 1)*(bstride) + ((bx) >> 1)])

class DaalaDecoder {
private:
  FILE *input;
  const char *path;

  ogg_page page;
  ogg_sync_state oy;
  ogg_stream_state os;

  daala_info di;
  daala_comment dc;
  daala_setup_info *dsi;
  daala_dec_ctx *dctx;

  bool readPage();
  bool readPacket(ogg_packet *packet);
  bool readHeaders();
public:
  od_img img;
  int frame;

  DaalaDecoder();
  ~DaalaDecoder();

  bool open(const char *path);
  bool step();
  void reset();

  int getWidth() const;
  int getHeight() const;

  int getFrameWidth() const;
  int getFrameHeight() const;
  int getRunningFrameCount() const;

  bool setBlockSizeBuffer(unsigned char *buf, size_t buf_sz);
  bool setBandFlagsBuffer(unsigned int *buf, size_t buf_sz);
};

bool DaalaDecoder::readPage() {
  while (ogg_sync_pageout(&oy, &page) != 1) {
    char *buffer = ogg_sync_buffer(&oy, 4096);
    if (buffer == NULL) {
      return false;
    }
    int bytes = fread(buffer, 1, 4096, input);
    // End of file
    if (bytes == 0) {
      return false;
    }
    if (ogg_sync_wrote(&oy, bytes) != 0) {
      return false;
    }
  }
  return true;
}

bool DaalaDecoder::readPacket(ogg_packet *packet) {
  while (ogg_stream_packetout(&os, packet) != 1) {
    if (!readPage()) {
      return false;
    }
    if (ogg_stream_pagein(&os, &page) != 0) {
      return false;
    }
  }
  return true;
}

bool DaalaDecoder::readHeaders() {
  bool done = false;
  while (!done && readPage()) {
    int serial = ogg_page_serialno(&page);
    if (ogg_page_bos(&page)) {
      if (ogg_stream_init(&os, serial) != 0) {
        return false;
      }
    }
    if (ogg_stream_pagein(&os, &page) != 0) {
      return false;
    }
    ogg_packet packet;
    while (!done && ogg_stream_packetpeek(&os, &packet) != 0) {
      int ret = daala_decode_header_in(&di, &dc, &dsi, &packet);
      if (ret < 0) {
        if (memcmp(packet.packet, "fishead", packet.bytes)) {
          fprintf(stderr, "Ogg Skeleton streams not supported\n");
        }
        return false;
      }
      if (ret == 0) {
        done = true;
        dctx = daala_decode_alloc(&di, dsi);
        if (dctx == NULL) {
          return false;
        }
      }
      if (!done) {
        if (ogg_stream_packetout(&os, &packet) != 1) {
          return false;
        }
      }
    }
  }
  return done;
}

DaalaDecoder::DaalaDecoder() : input(NULL), dsi(NULL), dctx(NULL) {
  daala_info_init(&di);
  daala_comment_init(&dc);
}

DaalaDecoder::~DaalaDecoder() {
  daala_info_clear(&di);
  daala_comment_clear(&dc);
  ogg_sync_clear(&oy);
  daala_setup_free(dsi);
  daala_decode_free(dctx);
}

bool DaalaDecoder::open(const char *path) {
  if (path == NULL) {
    return false;
  }
  ogg_sync_init(&oy);
  input = fopen(path,"rb");
  if (input == NULL) {
    return false;
  }
  this->path = path;
  frame = 0;
  return readHeaders();
}

bool DaalaDecoder::step() {
  //fprintf(stderr, "reading frame %i\n", frame);
  ogg_packet packet;
  if (readPacket(&packet)) {
    if (daala_decode_packet_in(dctx, &img, &packet) != 0) {
      return false;
    }
    frame++;
    return true;
  }
  return false;
}

// TODO: Move cleanup code from destructor to here
void DaalaDecoder::reset() {
}

int DaalaDecoder::getWidth() const {
  return di.pic_width;
}

int DaalaDecoder::getHeight() const {
  return di.pic_height;
}

int DaalaDecoder::getFrameWidth() const {
  return di.pic_width + (OD_BSIZE_MAX - 1) & ~(OD_BSIZE_MAX - 1);
}

int DaalaDecoder::getFrameHeight() const {
  return di.pic_height + (OD_BSIZE_MAX - 1) & ~(OD_BSIZE_MAX - 1);
}

bool DaalaDecoder::setBlockSizeBuffer(unsigned char *buf, size_t buf_sz) {
  if (dctx == NULL) {
    return false;
  }
  return daala_decode_ctl(dctx, OD_DECCTL_SET_BSIZE_BUFFER, buf, buf_sz) == 0;
}

bool DaalaDecoder::setBandFlagsBuffer(unsigned int *buf, size_t buf_sz) {
  if (dctx == NULL) {
    return false;
  }
  return daala_decode_ctl(dctx, OD_DECCTL_SET_FLAGS_BUFFER, buf, buf_sz) == 0;
}

int DaalaDecoder::getRunningFrameCount() const {
	return frame;
}

#define MIN_ZOOM (1)
#define MAX_ZOOM (4)

class TestPanel : public wxPanel {
  DECLARE_EVENT_TABLE()
private:
  DaalaDecoder dd;

  int zoom;
  unsigned char *pixels;

  unsigned char *bsize;
  int bstride;
  bool show_blocks;

  unsigned int *flags;
  int fstride;
  bool show_skip;
  bool show_noref;

  int getWidth() const;
  int getHeight() const;
  bool nextFrame();

  int getBand(int x, int y) const; 
public:
  TestPanel(wxWindow *parent);
  ~TestPanel();

  bool open(const char *path);
  void close();
  void render();

  int getZoom() const;
  bool setZoom(int zoom);

  void setShowBlocks(bool show_blocks);
  void setShowSkip(bool show_skip);
  void setShowNoRef(bool show_noref);


  void onKeyDown(wxKeyEvent &event);
  void onPaint(wxPaintEvent &event);
  void onIdle(wxIdleEvent &event);
  int getImageWidth() const;
  int getImageHeight() const;
};

class TestFrame : public wxFrame {
  DECLARE_EVENT_TABLE()
private:
  TestPanel *panel;
public:
  TestFrame();

  void onOpen(wxCommandEvent &event);
  void onClose(wxCommandEvent &event);
  void onQuit(wxCommandEvent &event);
  void onZoomIn(wxCommandEvent &event);
  void onZoomOut(wxCommandEvent &event);
  void onFilter(wxCommandEvent &event);
  void onAbout(wxCommandEvent &event);

  bool open(wxString path);
  void showFrameCount(int frame);
};

// Event table for TestPanel
BEGIN_EVENT_TABLE(TestPanel, wxPanel)
  EVT_KEY_DOWN(TestPanel::onKeyDown)
  EVT_PAINT(TestPanel::onPaint)
  //EVT_IDLE(TestPanel::onIdle)
END_EVENT_TABLE()

// Event table for TestFrame
enum {
  wxID_SHOW_BLOCKS = 6000,
  wxID_SHOW_SKIP,
  wxID_SHOW_NOREF
};

BEGIN_EVENT_TABLE(TestFrame, wxFrame)
  EVT_MENU(wxID_OPEN, TestFrame::onOpen)
  EVT_MENU(wxID_CLOSE, TestFrame::onClose)
  EVT_MENU(wxID_EXIT, TestFrame::onQuit)
  EVT_MENU(wxID_ZOOM_IN, TestFrame::onZoomIn)
  EVT_MENU(wxID_ZOOM_OUT, TestFrame::onZoomOut)
  EVT_MENU(wxID_SHOW_BLOCKS, TestFrame::onFilter)
  EVT_MENU(wxID_SHOW_SKIP, TestFrame::onFilter)
  EVT_MENU(wxID_SHOW_NOREF, TestFrame::onFilter)
  EVT_MENU(wxID_ABOUT, TestFrame::onAbout)
END_EVENT_TABLE()

// Enum for different aspects of Status display
enum {
  wxID_STATUS_FILE = 0,
  wxID_STATUS_FRAME_WIDTH,
  wxID_STATUS_FRAME_HEIGHT,
  wxID_STATUS_FRAME_COUNT
}

TestPanel::TestPanel(wxWindow *parent) : wxPanel(parent), pixels(NULL),
 zoom(0), bsize(NULL), show_blocks(false), flags(NULL), show_skip(false),
 show_noref(false) {
}

TestPanel::~TestPanel() {
  close();
}

bool TestPanel::open(const char *path) {
  if (!dd.open(path)) {
    return false;
  }
  if (!setZoom(MIN_ZOOM)) {
    return false;
  }
  int nhsb = dd.getFrameWidth() >> OD_LOG_BSIZE0 + OD_NBSIZES - 1;
  int nvsb = dd.getFrameHeight() >> OD_LOG_BSIZE0 + OD_NBSIZES - 1;
  bsize = (unsigned char *)malloc(sizeof(*bsize)*nhsb*4*nvsb*4);
  if (bsize == NULL) {
    close();
    return false;
  }
  bstride = nhsb*4;
  if (!dd.setBlockSizeBuffer(bsize, sizeof(*bsize)*nhsb*4*nvsb*4)) {
    close();
    return false;
  }
  flags = (unsigned int *)malloc(sizeof(*flags)*nhsb*8*nvsb*8);
  if (flags == NULL) {
    fprintf(stderr,"Could not allocate memory\n");
    close();
    return false;
  }
  fstride = nhsb*8;
  if (!dd.setBandFlagsBuffer(flags, sizeof(unsigned int)*nhsb*8*nvsb*8)) {
    fprintf(stderr,"Could not set flags buffer\n");
    close();
    return false;
  }
  if (!nextFrame()) {
    close();
    return false;
  }
  SetFocus();
  return true;
}

void TestPanel::close() {
  dd.reset();
  free(pixels);
  pixels = NULL;
  free(bsize);
  bsize = NULL;
  free(flags);
  flags = NULL;
}

int TestPanel::getWidth() const {
  return zoom*dd.getWidth();
}

int TestPanel::getHeight() const {
  return zoom*dd.getHeight();
}

int TestPanel::getBand(int x, int y) const {
  if (x == 0 && y == 0) return -1;
  if (x < 4 && y < 4) return 0;
  if (x < 8 && y < 2) return 1;
  if (x < 2 && y < 8) return 2;
  if (x < 8 && y < 8) return 3;
  if (x < 16 && y < 4) return 4;
  if (x < 4 && y < 16) return 5;
  if (x < 16 && y < 16) return 6;
  if (x < 32 && y < 8) return 7;
  if (x < 8 && y < 32) return 8;
  return 9;
}

int64_t block_edge_luma(int64_t yval) {
  return yval > 50 ? yval >> 1 : yval + 15;
}

void TestPanel::render() {
  od_img *img = &dd.img;
  /* Assume both chroma planes are decimated the same */
  int xdec = img->planes[1].xdec;
  int ydec = img->planes[1].ydec;
  int y_stride = img->planes[0].ystride;
  int cb_stride = img->planes[1].ystride;
  int cr_stride = img->planes[2].ystride;
  int p_stride = zoom*3*dd.getWidth();
  unsigned char *y_row = img->planes[0].data;
  unsigned char *cb_row = img->planes[1].data;
  unsigned char *cr_row = img->planes[2].data;
  unsigned char *p_row = pixels;
  for (int j = 0; j < dd.getHeight(); j++) {
    unsigned char *y = y_row;
    unsigned char *cb = cb_row;
    unsigned char *cr = cr_row;
    unsigned char *p = p_row;
    for (int i = 0; i < dd.getWidth(); i++) {
      int64_t yval;
      int64_t cbval;
      int64_t crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      yval = *y;
      cbval = *cb;
      crval = *cr;
      if (show_skip || show_noref) {
        unsigned char d = OD_BLOCK_SIZE4x4(bsize, bstride, i >> 2, j >> 2);
        int band = getBand(i & ((1 << d + 2) - 1), j & ((1 << d + 2) - 1));
        int bx = i & ~((1 << d + 2) - 1);
        int by = j & ~((1 << d + 2) - 1);
        unsigned int flag = flags[fstride*(by >> 2) + (bx >> 2)];
        cbval = 128;
        crval = 128;
        if (band >= 0) {
          /*R: U=84, V=255, B: U=255, V=107, G: U=43, V=21*/
          bool skip = (flag >> 2*band)&1;
          bool noref = (flag >> 2*band + 1)&1;
          if (skip && show_skip && noref && show_noref) {
            cbval = 43;
            crval = 21;
          }
          if ((!skip || !show_skip) && noref && show_noref) {
            cbval = 84;
            crval = 255;
          }
          if (skip && show_skip && (!noref || !show_noref)) {
            cbval = 255;
            crval = 107;
          }
        }
      }
      if (show_blocks) {
        unsigned char d = OD_BLOCK_SIZE4x4(bsize, bstride, i >> 2, j >> 2);
        int mask = (1 << d + OD_LOG_BSIZE0) - 1;
        if (!(i & mask) || !(j & mask)) {
          yval = block_edge_luma(yval);
          cbval = (cbval + 128) >> 1;
          crval = (crval + 128) >> 1;
        }
      }
      yval -= 16;
      cbval -= 128;
      crval -= 128;
      /*This is intentionally slow and very accurate.*/
      rval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 4490222169144LL*crval, 9745792000LL), 65535);
      gval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval - 534117096223LL*cbval - 1334761232047LL*crval,
       9745792000LL), 65535);
      bval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 5290866304968LL*cbval, 9745792000LL), 65535);
      unsigned char *px_row = p;
      for (int v = 0; v < zoom; v++) {
        unsigned char *px = px_row;
        for (int u = 0; u < zoom; u++) {
          *(px + 0) = (unsigned char)(rval >> 8);
          *(px + 1) = (unsigned char)(gval >> 8);
          *(px + 2) = (unsigned char)(bval >> 8);
          px += 3;
        }
        px_row += p_stride;
      }
      int dc = ((y - y_row) & 1) | (1 - xdec);
      y++;
      cb += dc;
      cr += dc;
      p += zoom*3;
    }
    int dc = -((j & 1) | (1 - ydec));
    y_row += y_stride;
    cb_row += dc & cb_stride;
    cr_row += dc & cr_stride;
    p_row += zoom*p_stride;
  }
}

int TestPanel::getZoom() const {
  return zoom;
}

bool TestPanel::setZoom(int z) {
  if (z <= MAX_ZOOM && z >= MIN_ZOOM && zoom != z) {
    unsigned char *p =
     (unsigned char *)malloc(sizeof(*p)*3*z*dd.getWidth()*z*dd.getHeight());
    if (p == NULL) {
      return false;
    }
    free(pixels);
    pixels = p;
    zoom = z;
    SetSize(getWidth(), getHeight());
    return true;
  }
  return false;
}

void TestPanel::setShowBlocks(bool show_blocks) {
  this->show_blocks = show_blocks;
}

void TestPanel::setShowSkip(bool show_skip) {
  this->show_skip = show_skip;
}

void TestPanel::setShowNoRef(bool show_noref) {
  this->show_noref = show_noref;
}

bool TestPanel::nextFrame() {
  if (dd.step()) {
    render();
    ((TestFrame* )GetParent())->showFrameCount(dd.getRunningFrameCount());
    return true;
  }
  return false;
}

int TestPanel::getImageWidth() const {
  return dd.getWidth();
}

int TestPanel::getImageHeight() const {
  return dd.getHeight();
}

void TestPanel::onKeyDown(wxKeyEvent &event) {
  switch (event.GetKeyCode()) {
    case '.' : {
      nextFrame();
      Refresh(false);
      break;
    }
  }
}

void TestPanel::onPaint(wxPaintEvent &) {
  wxBitmap bmp(wxImage(getWidth(), getHeight(), pixels, true));
  wxBufferedPaintDC dc(this, bmp);
}

void TestPanel::onIdle(wxIdleEvent &) {
  nextFrame();
  Refresh(false);
  /*wxMilliSleep(input.video_fps_n*1000/input.video_fps_n);*/
}

TestFrame::TestFrame() : wxFrame(NULL, wxID_ANY, _T("Daala Stream Analyzer"),
 wxDefaultPosition, wxDefaultSize, wxDEFAULT_FRAME_STYLE), panel(NULL) {
  wxMenuBar *mb = new wxMenuBar();
  wxMenu *fileMenu = new wxMenu();
  fileMenu->Append(wxID_OPEN, _T("&Open...\tAlt-O"), _T("Open daala file"));
  fileMenu->Append(wxID_CLOSE, _T("&Close\tAlt-C"), _T("Close daala file"));
  fileMenu->Enable(wxID_CLOSE, false);
  fileMenu->Append(wxID_EXIT, _T("E&xit\tAlt-X"), _T("Quit this program"));
  mb->Append(fileMenu, _T("&File"));
  wxMenu *viewMenu=new wxMenu();
  viewMenu->Append(wxID_ZOOM_IN, _T("Zoom-In\tCtrl-+"),
   _T("Double image size"));
  viewMenu->Append(wxID_ZOOM_OUT, _T("Zoom-Out\tCtrl--"),
   _T("Half image size"));
  viewMenu->AppendCheckItem(wxID_SHOW_BLOCKS, _T("&Blocks\tAlt-B"),
   _("Show block sizes"));
  viewMenu->AppendCheckItem(wxID_SHOW_SKIP, _T("&Skip\tAlt-S"),
   _("Show skip bands"));
  viewMenu->AppendCheckItem(wxID_SHOW_NOREF, _T("&No-Ref\tAlt-N"),
   _("Show no-ref bands"));
  mb->Append(viewMenu, _T("&View"));
  wxMenu *helpMenu=new wxMenu();
  helpMenu->Append(wxID_ABOUT, _T("&About...\tF1"), _T("Show about dialog"));
  mb->Append(helpMenu, _T("&Help"));
  SetMenuBar(mb);
  //File, Frame number, Width, Height
  CreateStatusBar(4);
  int status_widths[4] = {-1, -1, -1, -1};
  SetStatusWidths(4, status_widths);
  SetStatusText(_T("another day, another daala"));
}

void TestFrame::onOpen(wxCommandEvent& WXUNUSED(event)) {
  wxFileDialog openFileDialog(this, _T("Open file"), wxEmptyString,
   wxEmptyString, _T("Daala files (*.ogv)|*.ogv"),
   wxFD_OPEN | wxFD_FILE_MUST_EXIST);
  if (openFileDialog.ShowModal() != wxID_CANCEL) {
    open(openFileDialog.GetPath());
  }
}

void TestFrame::onClose(wxCommandEvent &WXUNUSED(event)) {
}

void TestFrame::onQuit(wxCommandEvent &WXUNUSED(event)) {
  Close(true);
}

void TestFrame::onZoomIn(wxCommandEvent &WXUNUSED(event)) {
  if (panel->setZoom(panel->getZoom() + 1)) {
    Fit();
    panel->render();
    panel->Refresh();
  }
}

void TestFrame::onZoomOut(wxCommandEvent &WXUNUSED(event)) {
  if (panel->setZoom(panel->getZoom() - 1)) {
    Fit();
    panel->render();
    panel->Refresh();
  }
}

void TestFrame::onFilter(wxCommandEvent &WXUNUSED(event)) {
  panel->setShowBlocks(GetMenuBar()->IsChecked(wxID_SHOW_BLOCKS));
  panel->setShowSkip(GetMenuBar()->IsChecked(wxID_SHOW_SKIP));
  panel->setShowNoRef(GetMenuBar()->IsChecked(wxID_SHOW_NOREF));
  panel->render();
  panel->Refresh(false);
}

void TestFrame::onAbout(wxCommandEvent& WXUNUSED(event)) {
  wxMessageBox(_T("asdf"), _T("basdf"), wxOK | wxICON_INFORMATION, this);
}

bool TestFrame::open(wxString path) {
  wxCharBuffer buffer = path.ToUTF8();
  const char *filename = buffer.data();
  panel = new TestPanel(this);
  if (panel->open(filename)) {
    Fit();
    SetStatusText(path, wxID_STATUS_FILE);
    SetStatusText(wxString::Format(wxT("W:%d"), panel->getImageWidth()), wxID_STATUS_FRAME_WIDTH);
    SetStatusText(wxString::Format(wxT("H:%d"), panel->getImageHeight()), wxID_STATUS_FRAME_HEIGHT);
    return true;
  }
  else {
    delete panel;
    panel = NULL;
    SetStatusText(_T("error loading file") + path, wxID_STATUS_FILE);
    return false;
  }
}

void TestFrame::showFrameCount(int frame) {
  SetStatusText(wxString::Format(wxT("Frame:%d"), frame-1), wxID_STATUS_FRAME_COUNT);
}

class TestApp : public wxApp {
private:
  TestFrame *frame;
public:
  bool OnInit();
};

bool TestApp::OnInit() {
  frame = new TestFrame();
  frame->Show();
  if (argc >= 2) {
    return frame->open(wxString(argv[1]));
  }
  return true;
}

IMPLEMENT_APP(TestApp)
