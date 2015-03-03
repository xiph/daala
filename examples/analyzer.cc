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

  bool setBlockSizeBuffer(unsigned char *buf, size_t buf_sz);
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

class TestPanel : public wxPanel {
  DECLARE_EVENT_TABLE()
private:
  DaalaDecoder dd;

  unsigned char *pixels;

  unsigned char *bsize;
  int bstride;
  bool show_blocks;

  int getWidth();
  int getHeight();
  bool nextFrame();
public:
  TestPanel(wxWindow *parent);
  ~TestPanel();

  bool open(const char *path);
  void close();
  void render();

  void setShowBlocks(bool show_blocks);

  void onKeyDown(wxKeyEvent &event);
  void onPaint(wxPaintEvent &event);
  void onIdle(wxIdleEvent &event);
};

BEGIN_EVENT_TABLE(TestPanel, wxPanel)
  EVT_KEY_DOWN(TestPanel::onKeyDown)
  EVT_PAINT(TestPanel::onPaint)
  //EVT_IDLE(TestPanel::onIdle)
END_EVENT_TABLE()

TestPanel::TestPanel(wxWindow *parent) : wxPanel(parent), pixels(NULL),
 bsize(NULL), show_blocks(false) {
}

TestPanel::~TestPanel() {
  close();
}

bool TestPanel::open(const char *path) {
  if (!dd.open(path)) {
    return false;
  }
  pixels =
   (unsigned char *)malloc(sizeof(*pixels)*3*dd.getWidth()*dd.getHeight());
  if (pixels == NULL) {
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
  if (!nextFrame()) {
    close();
    return false;
  }
  SetSize(getWidth(), getHeight());
  SetFocus();
  return true;
}

void TestPanel::close() {
  dd.reset();
  free(pixels);
  pixels = NULL;
  free(bsize);
  bsize = NULL;
}

int TestPanel::getWidth() {
  return dd.getWidth();
}

int TestPanel::getHeight() {
  return dd.getHeight();
}

void TestPanel::render() {
  od_img *img;
  int xdec;
  int ydec;
  unsigned char *y_row;
  unsigned char *cb_row;
  unsigned char *cr_row;
  unsigned char *y;
  unsigned char *cb;
  unsigned char *cr;
  int y_stride;
  int cb_stride;
  int cr_stride;
  unsigned char *p;
  int pitch;
  img = &dd.img;
  /* Assume both chroma planes are decimated the same */
  xdec = img->planes[1].xdec;
  ydec = img->planes[1].ydec;
  y_stride = img->planes[0].ystride;
  cb_stride = img->planes[1].ystride;
  cr_stride = img->planes[2].ystride;
  y_row = img->planes[0].data;
  cb_row = img->planes[1].data;
  cr_row = img->planes[2].data;
  p = pixels;
  pitch = 3*getWidth();
  for (int j = 0; j < getHeight(); j++) {
    int dc;
    y = y_row;
    cb = cb_row;
    cr = cr_row;
    for (int i = 0, k = 0; i < 3*getWidth(); i += 3, k++) {
      int64_t yval;
      int64_t cbval;
      int64_t crval;
      unsigned rval;
      unsigned gval;
      unsigned bval;
      yval = *y;
      if (show_blocks) {
        unsigned char d = OD_BLOCK_SIZE4x4(bsize, bstride, k >> 2, j >> 2);
        int mask = (1 << d + OD_LOG_BSIZE0) - 1;
        if (!(k & mask) || !(j & mask)) {
          yval >>= 1;
        }
      }
      yval -= 16;
      cbval = *cb - 128;
      crval = *cr - 128;
      /*This is intentionally slow and very accurate.*/
      rval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 4490222169144LL*crval, 9745792000LL), 65535);
      gval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval - 534117096223LL*cbval - 1334761232047LL*crval,
       9745792000LL), 65535);
      bval = OD_CLAMPI(0, (int32_t)OD_DIV_ROUND(
       2916394880000LL*yval + 5290866304968LL*cbval, 9745792000LL), 65535);
      *(p + pitch*j + i + 0) = (unsigned char)(rval >> 8);
      *(p + pitch*j + i + 1) = (unsigned char)(gval >> 8);
      *(p + pitch*j + i + 2) = (unsigned char)(bval >> 8);
      dc = ((y - y_row) & 1) | (1 - xdec);
      y++;
      cb += dc;
      cr += dc;
    }
    dc = -((j & 1) | (1 - ydec));
    y_row += y_stride;
    cb_row += dc & cb_stride;
    cr_row += dc & cr_stride;
  }
}

void TestPanel::setShowBlocks(bool show_blocks) {
  this->show_blocks = show_blocks;
}

bool TestPanel::nextFrame() {
  if (dd.step()) {
    render();
    return true;
  }
  return false;
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

class TestFrame : public wxFrame {
  DECLARE_EVENT_TABLE()
private:
  TestPanel *panel;
public:
  TestFrame();

  void onOpen(wxCommandEvent &event);
  void onClose(wxCommandEvent &event);
  void onQuit(wxCommandEvent &event);
  void onFilter(wxCommandEvent &event);
  void onAbout(wxCommandEvent &event);

  bool open(wxString path);
};

BEGIN_EVENT_TABLE(TestFrame, wxFrame)
  EVT_MENU(wxID_OPEN, TestFrame::onOpen)
  EVT_MENU(wxID_CLOSE, TestFrame::onClose)
  EVT_MENU(wxID_EXIT, TestFrame::onQuit)
  EVT_MENU(wxID_VIEW_DETAILS, TestFrame::onFilter)
  EVT_MENU(wxID_ABOUT, TestFrame::onAbout)
END_EVENT_TABLE()

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
  viewMenu->AppendCheckItem(wxID_VIEW_DETAILS, _T("&Blocks\tAlt-B"),
   _("Show blocks"));
  mb->Append(viewMenu, _T("&View"));
  wxMenu *helpMenu=new wxMenu();
  helpMenu->Append(wxID_ABOUT, _T("&About...\tF1"), _T("Show about dialog"));
  mb->Append(helpMenu, _T("&Help"));
  SetMenuBar(mb);
  CreateStatusBar(1);
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

void TestFrame::onFilter(wxCommandEvent &WXUNUSED(event)) {
  panel->setShowBlocks(GetMenuBar()->IsChecked(wxID_VIEW_DETAILS));
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
    wxBoxSizer *main = new wxBoxSizer(wxHORIZONTAL);
    main->Add(panel);
    SetSizerAndFit(main);
    SetStatusText(_T("loaded file: ") + path);
    return true;
  }
  else {
    delete panel;
    panel = NULL;
    SetStatusText(_T("error loading file") + path);
    return false;
  }
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
