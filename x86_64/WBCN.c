/* Created by Language version: 7.7.0 */
/* NOT VECTORIZED */
#define NRN_VECTORIZED 0
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib_ansi.h"
#undef PI
#define nil 0
#include "md1redef.h"
#include "section.h"
#include "nrniv_mf.h"
#include "md2redef.h"
 
#if METHOD3
extern int _method3;
#endif

#if !NRNGPU
#undef exp
#define exp hoc_Exp
extern double hoc_Exp(double);
#endif
 
#define nrn_init _nrn_init__WBCN
#define _nrn_initial _nrn_initial__WBCN
#define nrn_cur _nrn_cur__WBCN
#define _nrn_current _nrn_current__WBCN
#define nrn_jacob _nrn_jacob__WBCN
#define nrn_state _nrn_state__WBCN
#define _net_receive _net_receive__WBCN 
#define func func__WBCN 
#define states states__WBCN 
 
#define _threadargscomma_ /**/
#define _threadargsprotocomma_ /**/
#define _threadargs_ /**/
#define _threadargsproto_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define se _p[0]
#define gnabar _p[1]
#define gkbar _p[2]
#define gl _p[3]
#define el _p[4]
#define NNa _p[5]
#define NK _p[6]
#define PhiAs _p[7]
#define PhiBs _p[8]
#define il _p[9]
#define phiK _p[10]
#define phiNa _p[11]
#define n _p[12]
#define m _p[13]
#define h _p[14]
#define s _p[15]
#define qK _p[16]
#define pK _p[17]
#define qNa _p[18]
#define pNa _p[19]
#define ena _p[20]
#define ek _p[21]
#define ina _p[22]
#define ik _p[23]
#define am _p[24]
#define ah _p[25]
#define an _p[26]
#define as _p[27]
#define bm _p[28]
#define bh _p[29]
#define bn _p[30]
#define bs _p[31]
#define xiK _p[32]
#define xiNa _p[33]
#define Dn _p[34]
#define Dm _p[35]
#define Dh _p[36]
#define Ds _p[37]
#define DqK _p[38]
#define DpK _p[39]
#define DqNa _p[40]
#define DpNa _p[41]
#define _g _p[42]
#define _ion_ena	*_ppvar[0]._pval
#define _ion_ina	*_ppvar[1]._pval
#define _ion_dinadv	*_ppvar[2]._pval
#define _ion_ek	*_ppvar[3]._pval
#define _ion_ik	*_ppvar[4]._pval
#define _ion_dikdv	*_ppvar[5]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 
#if defined(__cplusplus)
extern "C" {
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static void _hoc_etam(void);
 static void _hoc_etah(void);
 static void _hoc_etan(void);
 static void _hoc_func(void);
 static int _mechtype;
extern void _nrn_cacheloop_reg(int, int);
extern void hoc_register_prop_size(int, int, int);
extern void hoc_register_limits(int, HocParmLimits*);
extern void hoc_register_units(int, HocParmUnits*);
extern void nrn_promote(Prop*, int, int);
extern Memb_func* memb_func;
 
#define NMODL_TEXT 1
#if NMODL_TEXT
static const char* nmodl_file_text;
static const char* nmodl_filename;
extern void hoc_reg_nmodl_text(int, const char*);
extern void hoc_reg_nmodl_filename(int, const char*);
#endif

 extern void _nrn_setdata_reg(int, void(*)(Prop*));
 static void _setdata(Prop* _prop) {
 _p = _prop->param; _ppvar = _prop->dparam;
 }
 static void _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range(int);
 _prop = hoc_getdata_range(_mechtype);
   _setdata(_prop);
 hoc_retpushx(1.);
}
 /* connect user functions to hoc names */
 static VoidFunc hoc_intfunc[] = {
 "setdata_WBCN", _hoc_setdata,
 "etam_WBCN", _hoc_etam,
 "etah_WBCN", _hoc_etah,
 "etan_WBCN", _hoc_etan,
 "func_WBCN", _hoc_func,
 0, 0
};
#define etam etam_WBCN
#define etah etah_WBCN
#define etan etan_WBCN
 extern double etam( double , double );
 extern double etah( double , double );
 extern double etan( double , double );
 /* declare global and static user variables */
#define TNa TNa_WBCN
 double TNa = 800;
#define TK TK_WBCN
 double TK = 400;
#define gamNa gamNa_WBCN
 double gamNa = 10;
#define gamK gamK_WBCN
 double gamK = 10;
#define tau tau_WBCN
 double tau = 1;
#define wNa2 wNa2_WBCN
 double wNa2 = 200;
#define wK2 wK2_WBCN
 double wK2 = 150;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 "gl_WBCN", 0, 1e+09,
 "gkbar_WBCN", 0, 1e+09,
 "gnabar_WBCN", 0, 1e+09,
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "tau_WBCN", "ms",
 "gnabar_WBCN", "S/cm2",
 "gkbar_WBCN", "S/cm2",
 "gl_WBCN", "S/cm2",
 "el_WBCN", "mV",
 "il_WBCN", "mA/cm2",
 0,0
};
 static double delta_t = 0.01;
 static double h0 = 0;
 static double m0 = 0;
 static double n0 = 0;
 static double pNa0 = 0;
 static double pK0 = 0;
 static double qNa0 = 0;
 static double qK0 = 0;
 static double s0 = 0;
 static double v = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "tau_WBCN", &tau_WBCN,
 "gamK_WBCN", &gamK_WBCN,
 "wK2_WBCN", &wK2_WBCN,
 "TK_WBCN", &TK_WBCN,
 "gamNa_WBCN", &gamNa_WBCN,
 "wNa2_WBCN", &wNa2_WBCN,
 "TNa_WBCN", &TNa_WBCN,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(Prop*);
static void  nrn_init(_NrnThread*, _Memb_list*, int);
static void nrn_state(_NrnThread*, _Memb_list*, int);
 static void nrn_cur(_NrnThread*, _Memb_list*, int);
static void  nrn_jacob(_NrnThread*, _Memb_list*, int);
 
static int _ode_count(int);
static void _ode_map(int, double**, double**, double*, Datum*, double*, int);
static void _ode_spec(_NrnThread*, _Memb_list*, int);
static void _ode_matsol(_NrnThread*, _Memb_list*, int);
 
#define _cvode_ieq _ppvar[6]._i
 static void _ode_matsol_instance1(_threadargsproto_);
 /* connect range variables in _p that hoc is supposed to know about */
 static const char *_mechanism[] = {
 "7.7.0",
"WBCN",
 "se_WBCN",
 "gnabar_WBCN",
 "gkbar_WBCN",
 "gl_WBCN",
 "el_WBCN",
 "NNa_WBCN",
 "NK_WBCN",
 "PhiAs_WBCN",
 "PhiBs_WBCN",
 0,
 "il_WBCN",
 "phiK_WBCN",
 "phiNa_WBCN",
 0,
 "n_WBCN",
 "m_WBCN",
 "h_WBCN",
 "s_WBCN",
 "qK_WBCN",
 "pK_WBCN",
 "qNa_WBCN",
 "pNa_WBCN",
 0,
 0};
 static Symbol* _na_sym;
 static Symbol* _k_sym;
 
extern Prop* need_memb(Symbol*);

static void nrn_alloc(Prop* _prop) {
	Prop *prop_ion;
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 43, _prop);
 	/*initialize range parameters*/
 	se = -1;
 	gnabar = 0.08;
 	gkbar = 0.036;
 	gl = 0.0003;
 	el = -54.3;
 	NNa = 6000;
 	NK = 1800;
 	PhiAs = 5e-05;
 	PhiBs = 0.00051;
 	_prop->param = _p;
 	_prop->param_size = 43;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 7, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_na_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[0]._pval = &prop_ion->param[0]; /* ena */
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ina */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dinadv */
 prop_ion = need_memb(_k_sym);
 nrn_promote(prop_ion, 0, 1);
 	_ppvar[3]._pval = &prop_ion->param[0]; /* ek */
 	_ppvar[4]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[5]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static void _initlists();
  /* some states have an absolute tolerance */
 static Symbol** _atollist;
 static HocStateTolerance _hoc_state_tol[] = {
 0,0
};
 static void _update_ion_pointer(Datum*);
 extern Symbol* hoc_lookup(const char*);
extern void _nrn_thread_reg(int, int, void(*)(Datum*));
extern void _nrn_thread_table_reg(int, void(*)(double*, Datum*, Datum*, _NrnThread*, int));
extern void hoc_register_tolerance(int, HocStateTolerance*, Symbol***);
extern void _cvode_abstol( Symbol**, double*, int);

 void _WBCN_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("na", -10000.);
 	ion_reg("k", -10000.);
 	_na_sym = hoc_lookup("na_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_setdata_reg(_mechtype, _setdata);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
 #if NMODL_TEXT
  hoc_reg_nmodl_text(_mechtype, nmodl_file_text);
  hoc_reg_nmodl_filename(_mechtype, nmodl_filename);
#endif
  hoc_register_prop_size(_mechtype, 43, 7);
  hoc_register_dparam_semantics(_mechtype, 0, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 1, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 2, "na_ion");
  hoc_register_dparam_semantics(_mechtype, 3, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 4, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 5, "k_ion");
  hoc_register_dparam_semantics(_mechtype, 6, "cvodeieq");
 	hoc_register_cvode(_mechtype, _ode_count, _ode_map, _ode_spec, _ode_matsol);
 	hoc_register_tolerance(_mechtype, _hoc_state_tol, &_atollist);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 WBCN /Users/giuliafranco/Desktop/Neuron_Simulations/mods/WBCN.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "WBCN.mod   squid sodium, potassium, and leak channels";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static void _modl_cleanup(){ _match_recurse=1;}
static int func(double);
 
static int _ode_spec1(_threadargsproto_);
/*static int _ode_matsol1(_threadargsproto_);*/
 static int _slist1[7], _dlist1[7];
 static int states(_threadargsproto_);
 
/*CVODE*/
 static int _ode_spec1 () {_reset=0;
 {
   double _lflag , _ln0 , _lh0 , _ls0 ;
 func ( _threadargscomma_ v ) ;
   _lflag = 0.0 ;
   _ln0 = n ;
   while ( _lflag  == 0.0 ) {
     Dn = 2.0 * ( an * ( 1.0 - n ) - bn * n + etan ( _threadargscomma_ an , bn ) ) ;
     if ( n < 0.0  || n > 1.0 ) {
       n = _ln0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
   _lflag = 0.0 ;
   _lh0 = h ;
   while ( _lflag  == 0.0 ) {
     Dh = 2.0 * ( ah * ( 1.0 - h ) - bh * h + etah ( _threadargscomma_ ah , bh ) ) ;
     if ( h < 0.0  || h > 1.0 ) {
       h = _lh0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
   _lflag = 0.0 ;
   _ls0 = s ;
   while ( _lflag  == 0.0 ) {
     Ds = as * ( 1.0 - s ) - bs * s ;
     if ( s < 0.0  || s > 1.0 ) {
       s = _ls0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
   DqK = pK / tau ;
   DpK = ( - gamK * pK / tau - wK2 * ( an * ( 1.0 - n ) + bn * n ) * qK ) + xiK ;
   DqNa = pNa / tau ;
   DpNa = ( - gamNa * pNa / tau - wNa2 * ( am * ( 1.0 - m ) + bm * m ) * qNa ) + xiNa ;
   }
 return _reset;
}
 static int _ode_matsol1 () {
 double _lflag , _ln0 , _lh0 , _ls0 ;
 func ( _threadargscomma_ v ) ;
 _lflag = 0.0 ;
 _ln0 = n ;
 while ( _lflag  == 0.0 ) {
   Dn = Dn  / (1. - dt*( ( 2.0 )*( ( ( an )*( ( ( - 1.0 ) ) ) - ( bn )*( 1.0 ) ) ) )) ;
   if ( n < 0.0  || n > 1.0 ) {
     n = _ln0 ;
     }
   else {
     _lflag = 1.0 ;
     }
   }
 _lflag = 0.0 ;
 _lh0 = h ;
 while ( _lflag  == 0.0 ) {
   Dh = Dh  / (1. - dt*( ( 2.0 )*( ( ( ah )*( ( ( - 1.0 ) ) ) - ( bh )*( 1.0 ) ) ) )) ;
   if ( h < 0.0  || h > 1.0 ) {
     h = _lh0 ;
     }
   else {
     _lflag = 1.0 ;
     }
   }
 _lflag = 0.0 ;
 _ls0 = s ;
 while ( _lflag  == 0.0 ) {
   Ds = Ds  / (1. - dt*( ( as )*( ( ( - 1.0 ) ) ) - ( bs )*( 1.0 ) )) ;
   if ( s < 0.0  || s > 1.0 ) {
     s = _ls0 ;
     }
   else {
     _lflag = 1.0 ;
     }
   }
 DqK = DqK  / (1. - dt*( 0.0 )) ;
 DpK = DpK  / (1. - dt*( ( ( ( - gamK )*( 1.0 ) ) / tau ) )) ;
 DqNa = DqNa  / (1. - dt*( 0.0 )) ;
 DpNa = DpNa  / (1. - dt*( ( ( ( - gamNa )*( 1.0 ) ) / tau ) )) ;
  return 0;
}
 /*END CVODE*/
 static int states () {_reset=0;
 {
   double _lflag , _ln0 , _lh0 , _ls0 ;
 func ( _threadargscomma_ v ) ;
   _lflag = 0.0 ;
   _ln0 = n ;
   while ( _lflag  == 0.0 ) {
      n = n + (1. - exp(dt*(( 2.0 )*( ( ( an )*( ( ( - 1.0 ) ) ) - ( bn )*( 1.0 ) ) ))))*(- ( ( 2.0 )*( ( ( an )*( ( 1.0 ) ) + etan ( _threadargscomma_ an , bn ) ) ) ) / ( ( 2.0 )*( ( ( an )*( ( ( - 1.0 ) ) ) - ( bn )*( 1.0 ) ) ) ) - n) ;
     if ( n < 0.0  || n > 1.0 ) {
       n = _ln0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
   _lflag = 0.0 ;
   _lh0 = h ;
   while ( _lflag  == 0.0 ) {
      h = h + (1. - exp(dt*(( 2.0 )*( ( ( ah )*( ( ( - 1.0 ) ) ) - ( bh )*( 1.0 ) ) ))))*(- ( ( 2.0 )*( ( ( ah )*( ( 1.0 ) ) + etah ( _threadargscomma_ ah , bh ) ) ) ) / ( ( 2.0 )*( ( ( ah )*( ( ( - 1.0 ) ) ) - ( bh )*( 1.0 ) ) ) ) - h) ;
     if ( h < 0.0  || h > 1.0 ) {
       h = _lh0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
   _lflag = 0.0 ;
   _ls0 = s ;
   while ( _lflag  == 0.0 ) {
      s = s + (1. - exp(dt*(( as )*( ( ( - 1.0 ) ) ) - ( bs )*( 1.0 ))))*(- ( ( as )*( ( 1.0 ) ) ) / ( ( as )*( ( ( - 1.0 ) ) ) - ( bs )*( 1.0 ) ) - s) ;
     if ( s < 0.0  || s > 1.0 ) {
       s = _ls0 ;
       }
     else {
       _lflag = 1.0 ;
       }
     }
    qK = qK - dt*(- ( ( pK ) / tau ) ) ;
    pK = pK + (1. - exp(dt*(( ( ( - gamK )*( 1.0 ) ) / tau ))))*(- ( ( ( - ( ( wK2 )*( ( ( an )*( ( 1.0 - n ) ) + ( bn )*( n ) ) ) )*( qK ) ) ) + xiK ) / ( ( ( ( - gamK )*( 1.0 ) ) / tau ) ) - pK) ;
    qNa = qNa - dt*(- ( ( pNa ) / tau ) ) ;
    pNa = pNa + (1. - exp(dt*(( ( ( - gamNa )*( 1.0 ) ) / tau ))))*(- ( ( ( - ( ( wNa2 )*( ( ( am )*( ( 1.0 - m ) ) + ( bm )*( m ) ) ) )*( qNa ) ) ) + xiNa ) / ( ( ( ( - gamNa )*( 1.0 ) ) / tau ) ) - pNa) ;
   }
  return 0;
}
 
double etan (  double _lan , double _lbn ) {
   double _letan;
  _letan = normrand ( 0.0 , sqrt ( _lan * ( 1.0 - n ) + _lbn * n ) * dt / sqrt ( 4.0 * NK ) ) ;
    
return _letan;
 }
 
static void _hoc_etan(void) {
  double _r;
   _r =  etan (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double etam (  double _lam , double _lbm ) {
   double _letam;
  _letam = normrand ( 0.0 , sqrt ( _lam * ( 1.0 - m ) + _lbm * m ) * dt / sqrt ( 3.0 * NNa ) ) ;
    
return _letam;
 }
 
static void _hoc_etam(void) {
  double _r;
   _r =  etam (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
double etah (  double _lah , double _lbh ) {
   double _letah;
  _letah = normrand ( 0.0 , sqrt ( _lah * ( 1.0 - h ) + _lbh * h ) * dt / sqrt ( NNa ) ) ;
    
return _letah;
 }
 
static void _hoc_etah(void) {
  double _r;
   _r =  etah (  *getarg(1) , *getarg(2) );
 hoc_retpushx(_r);
}
 
static int  func (  double _lv ) {
   double _lq10 ;
  as = PhiAs * exp ( - ( _lv + 85. ) / 30. ) ;
   bs = PhiBs / ( exp ( - 0.3 * ( _lv + 17. ) ) + 1.0 ) ;
   am = 0.1 * ( _lv + 40.0 ) / ( 1.0 - exp ( - ( _lv + 40.0 ) / 10.0 ) ) ;
   bm = 4.0 * exp ( - ( _lv + 65.0 ) / 18.0 ) ;
   ah = 0.07 * exp ( - ( _lv + 65.0 ) / 20.0 ) ;
   bh = 1.0 / ( 1.0 + exp ( - ( _lv + 35.0 ) / 10.0 ) ) ;
   an = 0.01 * ( _lv + 55.0 ) / ( 1.0 - exp ( - ( _lv + 55.0 ) / 10.0 ) ) ;
   bn = 0.125 * exp ( - ( _lv + 65.0 ) / 80.0 ) ;
   xiK = normrand ( 0.0 , sqrt ( gamK * TK * ( an * ( 1.0 - n ) + bn * n ) ) * dt ) ;
   xiNa = normrand ( 0.0 , sqrt ( gamNa * TNa * ( am * ( 1.0 - m ) + bm * m ) ) * dt ) ;
   phiK = sqrt ( ( pow( n , 4.0 ) * ( 1.0 - pow( n , 4.0 ) ) ) / NK ) * qK ;
   phiNa = sqrt ( ( pow( m , 3.0 ) * ( 1.0 - pow( m , 3.0 ) ) ) / NNa ) * h * qNa ;
     return 0; }
 
static void _hoc_func(void) {
  double _r;
   _r = 1.;
 func (  *getarg(1) );
 hoc_retpushx(_r);
}
 
static int _ode_count(int _type){ return 7;}
 
static void _ode_spec(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
     _ode_spec1 ();
   }}
 
static void _ode_map(int _ieq, double** _pv, double** _pvdot, double* _pp, Datum* _ppd, double* _atol, int _type) { 
 	int _i; _p = _pp; _ppvar = _ppd;
	_cvode_ieq = _ieq;
	for (_i=0; _i < 7; ++_i) {
		_pv[_i] = _pp + _slist1[_i];  _pvdot[_i] = _pp + _dlist1[_i];
		_cvode_abstol(_atollist, _atol, _i);
	}
 }
 
static void _ode_matsol_instance1(_threadargsproto_) {
 _ode_matsol1 ();
 }
 
static void _ode_matsol(_NrnThread* _nt, _Memb_list* _ml, int _type) {
   Datum* _thread;
   Node* _nd; double _v; int _iml, _cntml;
  _cntml = _ml->_nodecount;
  _thread = _ml->_thread;
  for (_iml = 0; _iml < _cntml; ++_iml) {
    _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
    _nd = _ml->_nodelist[_iml];
    v = NODEV(_nd);
  ena = _ion_ena;
  ek = _ion_ek;
 _ode_matsol_instance1(_threadargs_);
 }}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_na_sym, _ppvar, 0, 0);
   nrn_update_ion_pointer(_na_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_na_sym, _ppvar, 2, 4);
   nrn_update_ion_pointer(_k_sym, _ppvar, 3, 0);
   nrn_update_ion_pointer(_k_sym, _ppvar, 4, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 5, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  h = h0;
  m = m0;
  n = n0;
  pNa = pNa0;
  pK = pK0;
  qNa = qNa0;
  qK = qK0;
  s = s0;
 {
   func ( _threadargscomma_ v ) ;
   if ( se > 0.0 ) {
     set_seed ( se ) ;
     }
   n = an / ( an + bn ) ;
   m = am / ( am + bm ) ;
   h = ah / ( ah + bh ) ;
   s = as / ( as + bs ) ;
   qK = 0.0 ;
   pK = 0.0 ;
   qNa = 0.0 ;
   pNa = 0.0 ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  ena = _ion_ena;
  ek = _ion_ek;
 initmodel();
  }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   m = am / ( am + bm ) ;
   ina = gnabar * ( pow( m , 3.0 ) * h * s + phiNa ) * ( v - ena ) ;
   ik = gkbar * ( pow( n , 4.0 ) + phiK ) * ( v - ek ) ;
   il = gl * ( v - el ) ;
   }
 _current += ina;
 _current += ik;
 _current += il;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  ena = _ion_ena;
  ek = _ion_ek;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
 double _dina;
  _dina = ina;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dinadv += (_dina - ina)/.001 ;
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ina += ina ;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v = 0.0; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v=_v;
{
  ena = _ion_ena;
  ek = _ion_ek;
 { error =  states();
 if(error){fprintf(stderr,"at line 78 in file WBCN.mod:\n	SOLVE states METHOD cnexp\n"); nrn_complain(_p); abort_run(error);}
 }  }}

}

static void terminal(){}

static void _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
 _slist1[0] = &(n) - _p;  _dlist1[0] = &(Dn) - _p;
 _slist1[1] = &(h) - _p;  _dlist1[1] = &(Dh) - _p;
 _slist1[2] = &(s) - _p;  _dlist1[2] = &(Ds) - _p;
 _slist1[3] = &(qK) - _p;  _dlist1[3] = &(DqK) - _p;
 _slist1[4] = &(pK) - _p;  _dlist1[4] = &(DpK) - _p;
 _slist1[5] = &(qNa) - _p;  _dlist1[5] = &(DqNa) - _p;
 _slist1[6] = &(pNa) - _p;  _dlist1[6] = &(DpNa) - _p;
_first = 0;
}

#if NMODL_TEXT
static const char* nmodl_filename = "/Users/giuliafranco/Desktop/Neuron_Simulations/mods/WBCN.mod";
static const char* nmodl_file_text = 
  "TITLE WBCN.mod   squid sodium, potassium, and leak channels\n"
  " \n"
  "COMMENT\n"
  "Stochastic Wang Buzsaki equations with colored noise (WBCN)\n"
  "Equations as in Guler (2013) Neural Comp. 25:46-74\n"
  "Implemented for Pezo, Soudry and Orio (2014) Front Comp Neurosci\n"
  "\n"
  "ENDCOMMENT\n"
  " \n"
  "UNITS {\n"
  "	(mA) = (milliamp)\n"
  "	(mV) = (millivolt)\n"
  "	(S) = (siemens)\n"
  "}\n"
  " \n"
  "NEURON {\n"
  "	SUFFIX WBCN\n"
  "	USEION na READ ena WRITE ina\n"
  "	USEION k READ ek WRITE ik\n"
  "	NONSPECIFIC_CURRENT il\n"
  "	RANGE gnabar, gkbar, gl, el, NNa, NK, se, m, h, n, s, phiNa, phiK, PhiAs, PhiBs\n"
  "}\n"
  " \n"
  "PARAMETER {\n"
  "    se = -1\n"
  "	gnabar = .08 (S/cm2)	<0,1e9>\n"
  "	gkbar = .036 (S/cm2)	<0,1e9>\n"
  "	gl = .0003 (S/cm2)	<0,1e9>\n"
  "	el = -54.3 (mV)\n"
  "	NNa = 6000\n"
  "	NK = 1800 \n"
  "	tau = 1 (ms)\n"
  "	gamK = 10\n"
  "	wK2 = 150\n"
  "	TK = 400\n"
  "	gamNa = 10\n"
  "	wNa2 = 200\n"
  "	TNa = 800\n"
  "	PhiAs = 0.00005\n"
  "	PhiBs = 0.00051\n"
  "}\n"
  " \n"
  "ASSIGNED {\n"
  "	v (mV)\n"
  "	celsius (degC)\n"
  "	ena (mV)\n"
  "	ek (mV)\n"
  "	dt (ms)\n"
  "	ina (mA/cm2)\n"
  "	ik (mA/cm2)\n"
  "	il (mA/cm2)\n"
  "	am	(/ms)\n"
  "	ah	(/ms)\n"
  "	an	(/ms)\n"
  "	as	(/ms)\n"
  "	bm	(/ms)\n"
  "	bh	(/ms)\n"
  "	bn	(/ms)\n"
  "	bs	(/ms)\n"
  "	phiK\n"
  "	phiNa\n"
  "	xiK	(/ms)\n"
  "	xiNa	(/ms)\n"
  "}\n"
  " \n"
  "STATE {	\n"
  "	n\n"
  "	m\n"
  "	h\n"
  "	s\n"
  "	qK\n"
  "	pK\n"
  "	qNa\n"
  "	pNa\n"
  "}\n"
  "\n"
  "BREAKPOINT {\n"
  "	SOLVE states METHOD cnexp\n"
  "    m=am/(am+bm)\n"
  "	ina = gnabar*(m^3*h*s+phiNa)*(v - ena)\n"
  "	ik = gkbar*(n^4+phiK)*(v - ek)\n"
  "	il = gl*(v - el)\n"
  "}\n"
  " \n"
  "INITIAL {\n"
  "	func(v)\n"
  "    if (se>0) {set_seed(se)} \n"
  "    n=an/(an+bn)\n"
  "	m=am/(am+bm)\n"
  "	h=ah/(ah+bh)\n"
  "	s=as/(as+bs)\n"
  "	qK=0\n"
  "	pK=0\n"
  "	qNa=0\n"
  "	pNa=0\n"
  "}\n"
  "\n"
  "DERIVATIVE states {  \n"
  "	LOCAL flag, n0, h0, s0\n"
  "	func(v)\n"
  "\n"
  "	flag=0\n"
  "	n0=n\n"
  "    : If the variable leaves [0,1] then random numbers are drawn again\n"
  "	while (flag==0) {\n"
  "		n'=2*(an*(1-n)-bn*n+etan(an,bn))\n"
  "		if (n<0 || n>1) {\n"
  "			n=n0\n"
  "		} else {flag=1}\n"
  "	}\n"
  "	\n"
  "	flag=0\n"
  "	h0=h\n"
  "	while (flag==0) {\n"
  "		h'=2*(ah*(1-h)-bh*h+etah(ah,bh))\n"
  "		if (h<0 || h>1) {\n"
  "			h=h0\n"
  "		} else {flag=1}\n"
  "	}\n"
  "	flag=0\n"
  "	s0=s\n"
  "	while (flag==0) {\n"
  "		s'=as*(1-s)-bs*s\n"
  "		if (s<0 || s>1) {\n"
  "			s=s0\n"
  "		} else {flag=1}\n"
  "	}\n"
  "\n"
  "	qK'=pK/tau\n"
  "	pK'=(-gamK*pK/tau-wK2*(an*(1-n)+bn*n)*qK)+xiK\n"
  "	\n"
  "	qNa'=pNa/tau\n"
  "	pNa'=(-gamNa*pNa/tau-wNa2*(am*(1-m)+bm*m)*qNa)+xiNa\n"
  "}\n"
  "\n"
  "FUNCTION etan (an (/ms), bn (/ms)) (/ms) {\n"
  "	UNITSOFF\n"
  "	etan = normrand(0,sqrt(an*(1-n)+bn*n)*dt/sqrt(4*NK))\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION etam (am (/ms), bm (/ms)) (/ms) {\n"
  "	UNITSOFF\n"
  "	etam = normrand(0,sqrt(am*(1-m)+bm*m)*dt/sqrt(3*NNa))\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "FUNCTION etah (ah (/ms), bh (/ms)) (/ms) {\n"
  "	UNITSOFF\n"
  "	etah = normrand(0,sqrt(ah*(1-h)+bh*h)*dt/sqrt(NNa))\n"
  "	UNITSON\n"
  "}\n"
  "\n"
  "PROCEDURE func(v(mV)) {  :Computes rate and other constants at current v.\n"
  "	LOCAL q10\n"
  "	UNITSOFF\n"
  "	as=PhiAs*exp(-(v+85.)/30.)\n"
  "	bs=PhiBs/(exp(-0.3*(v+17.))+1)\n"
  "	am = 0.1*(v+40)/(1-exp(-(v+40)/10))\n"
  "	bm = 4*exp(-(v+65)/18)\n"
  "	ah = 0.07*exp(-(v+65)/20) \n"
  "	bh = 1/(1+exp(-(v+35)/10))\n"
  "	an = 0.01*(v+55)/(1-exp(-(v+55)/10))\n"
  "	bn = 0.125*exp(-(v+65)/80)\n"
  "\n"
  "	xiK=normrand(0,sqrt(gamK*TK*(an*(1-n)+bn*n))*dt)\n"
  "	xiNa=normrand(0,sqrt(gamNa*TNa*(am*(1-m)+bm*m))*dt)\n"
  "\n"
  "	phiK=sqrt((n^4*(1-n^4))/NK)*qK\n"
  "	phiNa=sqrt((m^3*(1-m^3))/NNa)*h*qNa\n"
  "	UNITSON \n"
  "}\n"
  ;
#endif
