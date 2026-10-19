"""Minimal directed MPFR bindings via ctypes (libmpfr.so.6 + libgmp.so.10)."""
from __future__ import annotations
import ctypes
from ctypes import c_int, c_long, c_double, c_char_p, c_void_p, POINTER, Structure, byref

gmp = ctypes.CDLL("libgmp.so.10")
mpfr = ctypes.CDLL("libmpfr.so.6")

MPFR_RNDN, MPFR_RNDZ, MPFR_RNDU, MPFR_RNDD = 0, 1, 2, 3  # nearest, zero, up, down

class __mpfr_struct(Structure):
    _fields_ = [
        ("_mpfr_prec", c_long),
        ("_mpfr_sign", c_int),
        ("_mpfr_exp", c_long),
        ("_mpfr_d", c_void_p),
    ]

mpfr_t = __mpfr_struct
mpfr.mpfr_init2.argtypes = [POINTER(mpfr_t), c_long]
mpfr.mpfr_clear.argtypes = [POINTER(mpfr_t)]
mpfr.mpfr_set_si.argtypes = [POINTER(mpfr_t), c_long, c_int]
mpfr.mpfr_set_ui.argtypes = [POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_set_d.argtypes = [POINTER(mpfr_t), c_double, c_int]
mpfr.mpfr_set_str.argtypes = [POINTER(mpfr_t), c_char_p, c_int, c_int]
mpfr.mpfr_set.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_add.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_sub.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_mul.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_div.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_mul_si.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_long, c_int]
mpfr.mpfr_mul_ui.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_div_ui.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_add_ui.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_sub_ui.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_ui_div.argtypes = [POINTER(mpfr_t), ctypes.c_ulong, POINTER(mpfr_t), c_int]
mpfr.mpfr_ui_pow_ui.argtypes = [POINTER(mpfr_t), ctypes.c_ulong, ctypes.c_ulong, c_int]
mpfr.mpfr_pow_ui.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), ctypes.c_ulong, c_int]
mpfr.mpfr_pow_si.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_long, c_int]
mpfr.mpfr_log.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_sqrt.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_neg.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t), c_int]
mpfr.mpfr_get_d.argtypes = [POINTER(mpfr_t), c_int]
mpfr.mpfr_get_d.restype = c_double
mpfr.mpfr_cmp.argtypes = [POINTER(mpfr_t), POINTER(mpfr_t)]
mpfr.mpfr_cmp.restype = c_int
mpfr.mpfr_snprintf = mpfr.mpfr_snprintf
# use mpfr_get_str for output
mpfr.mpfr_get_str.argtypes = [c_char_p, POINTER(c_long), c_int, ctypes.c_size_t, POINTER(mpfr_t), c_int]
mpfr.mpfr_get_str.restype = c_void_p
gmp.__gmp_free_func = None
# mpfr_free_str
mpfr.mpfr_free_str.argtypes = [c_void_p]

class Mpfr:
    __slots__ = ("_v",)
    def __init__(self, prec=256):
        self._v = mpfr_t()
        mpfr.mpfr_init2(byref(self._v), prec)
    def __del__(self):
        try:
            mpfr.mpfr_clear(byref(self._v))
        except Exception:
            pass
    def set_si(self, n, rnd=MPFR_RNDN):
        mpfr.mpfr_set_si(byref(self._v), n, rnd); return self
    def set_ui(self, n, rnd=MPFR_RNDN):
        mpfr.mpfr_set_ui(byref(self._v), n, rnd); return self
    def set_d(self, x, rnd=MPFR_RNDN):
        mpfr.mpfr_set_d(byref(self._v), float(x), rnd); return self
    def set_str(self, s, rnd=MPFR_RNDN, base=10):
        rc = mpfr.mpfr_set_str(byref(self._v), s.encode(), base, rnd)
        if rc != 0:
            raise ValueError(f"mpfr_set_str failed on {s!r}")
        return self
    def set(self, other, rnd=MPFR_RNDN):
        mpfr.mpfr_set(byref(self._v), byref(other._v), rnd); return self
    def add(self, a, b, rnd):
        mpfr.mpfr_add(byref(self._v), byref(a._v), byref(b._v), rnd); return self
    def sub(self, a, b, rnd):
        mpfr.mpfr_sub(byref(self._v), byref(a._v), byref(b._v), rnd); return self
    def mul(self, a, b, rnd):
        mpfr.mpfr_mul(byref(self._v), byref(a._v), byref(b._v), rnd); return self
    def div(self, a, b, rnd):
        mpfr.mpfr_div(byref(self._v), byref(a._v), byref(b._v), rnd); return self
    def mul_ui(self, a, n, rnd):
        mpfr.mpfr_mul_ui(byref(self._v), byref(a._v), n, rnd); return self
    def mul_si(self, a, n, rnd):
        mpfr.mpfr_mul_si(byref(self._v), byref(a._v), n, rnd); return self
    def div_ui(self, a, n, rnd):
        mpfr.mpfr_div_ui(byref(self._v), byref(a._v), n, rnd); return self
    def add_ui(self, a, n, rnd):
        mpfr.mpfr_add_ui(byref(self._v), byref(a._v), n, rnd); return self
    def sub_ui(self, a, n, rnd):
        mpfr.mpfr_sub_ui(byref(self._v), byref(a._v), n, rnd); return self
    def ui_div(self, n, a, rnd):
        mpfr.mpfr_ui_div(byref(self._v), n, byref(a._v), rnd); return self
    def pow_ui(self, a, n, rnd):
        mpfr.mpfr_pow_ui(byref(self._v), byref(a._v), n, rnd); return self
    def log(self, a, rnd):
        mpfr.mpfr_log(byref(self._v), byref(a._v), rnd); return self
    def sqrt(self, a, rnd):
        mpfr.mpfr_sqrt(byref(self._v), byref(a._v), rnd); return self
    def neg(self, a, rnd):
        mpfr.mpfr_neg(byref(self._v), byref(a._v), rnd); return self
    def get_d(self, rnd=MPFR_RNDN):
        return mpfr.mpfr_get_d(byref(self._v), rnd)
    def to_str(self, digits=40, rnd=MPFR_RNDN):
        exp = c_long()
        p = mpfr.mpfr_get_str(None, byref(exp), 10, digits, byref(self._v), rnd)
        if not p:
            raise RuntimeError("mpfr_get_str failed")
        s = ctypes.string_at(p).decode()
        mpfr.mpfr_free_str(p)
        # s is mantissa digits; exp is 1-based exponent for leading digit
        if s.startswith("-"):
            sign, digits_s = "-", s[1:]
        else:
            sign, digits_s = "", s
        e = exp.value
        if e <= 0:
            return f"{sign}0.{'0'*(-e)}{digits_s}"
        if e >= len(digits_s):
            return f"{sign}{digits_s}{'0'*(e-len(digits_s))}"
        return f"{sign}{digits_s[:e]}.{digits_s[e:]}"
