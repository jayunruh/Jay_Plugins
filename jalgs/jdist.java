/*******************************************************************************
 * Copyright (c) 2013 Jay Unruh, Stowers Institute for Medical Research.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v2.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/old-licenses/gpl-2.0.html
 ******************************************************************************/

package jalgs;

public class jdist{
	// here we have utility methods for calculating limiting statistical values
	// most of the routines are adapted from numerical recipes

	public double ibeta(double a,double b,double x){
		// Returns the incomplete beta function Ix(a, b). Adapted from numerical
		// recipes.
		double beta=0.0;
		if(x<0.0||x>1.0){
			return Double.NaN;
		}
		if(x!=0.0||x!=1.0){
			beta=Math.exp(sf.gammaln(a+b)-sf.gammaln(a)-sf.gammaln(b)+a*Math.log(x)+b*Math.log(1.0-x));
		}
		if(x<(a+1.0)/(a+b+2.0)){
			return beta*betacf(a,b,x)/a;
		}else{
			return 1.0-beta*betacf(b,a,1.0-x)/b;
		}
	}

	private double betacf(double a,double b,double x){
		// adapted from numerical recipes
		int MAXIT=100;
		double EPS=3.0e-7;
		double FPMIN=1.0e-30;
		double qab=a+b;
		double qap=a+1.0;
		double qam=a-1.0;
		double c=1.0;
		double d=1.0-qab*x/qap;
		if(Math.abs(d)<FPMIN){
			d=FPMIN;
		}
		d=1.0/d;
		double h=d;
		int m=0;
		while(m<MAXIT){
			int m2=2*(m+1);
			double aa=(m+1)*(b-m-1)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if(Math.abs(d)<FPMIN){
				d=FPMIN;
			}
			c=1.0+aa/c;
			if(Math.abs(c)<FPMIN){
				c=FPMIN;
			}
			d=1.0/d;
			h*=d*c;
			aa=-(a+m+1)*(qab+m+1)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if(Math.abs(d)<FPMIN)
				d=FPMIN;
			c=1.0+aa/c;
			if(Math.abs(c)<FPMIN)
				c=FPMIN;
			d=1.0/d;
			double del=d*c;
			h*=del;
			if(Math.abs(del-1.0)<EPS){
				return h;
			}
			m++;
		}
		return Double.NaN;
	}

	public double FCumProb(double numdof,double dendof,double F){
		double xval=numdof*F/(numdof*F+dendof);
		return ibeta(0.5*numdof,0.5*dendof,xval);
	}

	public double FLimit(double numdof,double dendof,double prob){
		double df=0.001;
		double F=1.0-df;
		double cumprob=0.0;
		double lastcumprob=0.0;
		do{
			F+=df;
			lastcumprob=cumprob;
			cumprob=FCumProb(numdof,dendof,F);
		}while(cumprob<prob);
		double fraction=(prob-lastcumprob)/(cumprob-lastcumprob);
		return fraction*(F-df)+(1.0-fraction)*F;
	}

	public double tCumProbTwoSampEqVar(double dof1,double dof2,double mean1,double mean2,double sterr1,double sterr2,double nulldiff,boolean twotailed){
		double n1=dof1+1.0;
		double n2=dof2+1.0;
		// assume that variances are equal and that dof=n-1;
		double S2pooled=(dof1*n1*sterr1*sterr1+dof2*n2*sterr2*sterr2)/(dof1+dof2);
		double sterrdiff=Math.sqrt(S2pooled*(1.0/n1+1.0/n2));
		double t=(Math.abs(mean1-mean2)-nulldiff)/sterrdiff;
		return tCumProb(dof1+dof2,t,twotailed);
	}

	public double tCumProbTwoSampNVar(double dof1,double dof2,double mean1,double mean2,double sterr1,double sterr2,double nulldiff,boolean twotailed){
		double n1=dof1+1.0;
		double n2=dof2+1.0;
		// assume that variances are not equal and that dof=n-1;
		double sterrdiff=Math.sqrt(sterr1*sterr1+sterr2*sterr2);
		double C=sterr1*sterr1/(sterr1*sterr1+sterr2*sterr2);
		double dof=dof1*dof2/(dof2*C*C+dof1*(1.0-C)*(1.0-C));
		double t=(Math.abs(mean1-mean2)-nulldiff)/sterrdiff;
		return tCumProb(dof,t,twotailed);
	}

	public double tCumProbOneSamp(double dof,double mean,double value,double sterr,boolean twotailed){
		// t is the number of standard errors the value is away from the mean
		// the standard error is given by stdev/sqrt(n);
		double t=(value-mean)/sterr;
		return tCumProb(dof,t,twotailed);
	}
	
	public double tLim(double dof,double prob,boolean twotailed){
		//this gets the critical t value for a specific probability
		double dt=0.001;
		double t=0.0-dt;
		double cumprob=0.0;
		double lastcumprob=0.0;
		do{
			t+=dt;
			lastcumprob=cumprob;
			cumprob=tCumProb(dof,t,twotailed);
		}while(cumprob>prob);
		double fraction=(prob-cumprob)/(lastcumprob-cumprob);
		return t-fraction*dt;
	}

	public double tCumProb(double dof,double t,boolean twotailed){
		//this gets the probabity for a specific t value
		if(!twotailed){
			return 0.5*ibeta(0.5*dof,0.5,dof/(dof+t*t));
		}else{
			return ibeta(0.5*dof,0.5,dof/(dof+t*t));
		}
	}
	
	/**
	 * this code is copied from jdistlib under the gnu gpl license
	 * Perform Hartigan's dip test, assuming the minimum test statistics D is zero.
	 * @param x MUST BE SORTED in order to output the right result. This routine will NOT check for order!
	 * @return an array of four elements: The first is the test statistic, the second is the p-value, followed by indices for which there are a dip. If there is no dip, the indices will be set to -1.
	 */
	public static final double[] diptest_presorted(double[] x) {
		double dip = 0, dip_l, dip_u, dipnew;
		int n = x.length, low = 1, high = n, l_gcm = 0, l_lcm = 0;
		int mnj, mnmnj, mjk, mjmjk, ig, ih, iv, ix, i;
		int[] gcm = new int[n+1], lcm = new int[n+1], mn = new int[n+1], mj = new int[n+1];
		if (n < 2 || x[0] == x[n-1]) return new double[] {0, 1, -1, -1};
		//* Establish the indices  mn[0..n-1]  over which combination is necessary for the convex MINORANT (GCM) fit.
		mn[1] = 1;
		for (int j = 2; j <= n; ++j) {
			mn[j] = j - 1;
			while(true) {
				mnj = mn[j];
				mnmnj = mn[mnj];
				if (mnj == 1 || ( x[j-1]  - x[mnj-1]) * (mnj - mnmnj) < (x[mnj-1] - x[mnmnj-1]) * (j - mnj)) break;
				mn[j] = mnmnj;
			}
		}
		// Establish the indices   mj[0..n-1]  over which combination is necessary for the concave MAJORANT (LCM) fit.
		mj[n] = n;
		for (int k = n - 1; k >= 1; k--) {
			mj[k] = k + 1;
			while(true) {
				mjk = mj[k];
				mjmjk = mj[mjk];
				if (mjk == n || ( x[k-1]  - x[mjk-1]) * (mjk - mjmjk) < (x[mjk-1] - x[mjmjk-1]) * (k - mjk)) break;
				mj[k] = mjmjk;
			}
		}

		/* ----------------------- Start the cycling. ------------------------------- */
		//LOOP_Start:
		while (true) {

			/* Collect the change points for the GCM from HIGH to LOW. */
			gcm[1] = high;
			for(i = 1; gcm[i] > low; i++)
				gcm[i+1] = mn[gcm[i]];
			ig = l_gcm = i; // l_gcm == relevant_length(GCM)
			ix = ig - 1; //  ix, ig  are counters for the convex minorant.

			/* Collect the change points for the LCM from LOW to HIGH. */
			lcm[1] = low;
			for(i = 1; lcm[i] < high; i++)
				lcm[i+1] = mj[lcm[i]];
			ih = l_lcm = i; // l_lcm == relevant_length(LCM)
			iv = 2; //  iv, ih  are counters for the concave majorant.

			//	Find the largest distance greater than 'DIP' between the GCM and the LCM from LOW to HIGH.

			// FIXME: <Rconfig.h>  should provide LDOUBLE or something like it
			/* long */ double d = 0.;// <<-- see if this makes 32-bit/64-bit difference go..
			if (l_gcm != 2 || l_lcm != 2) {
				//if(*debug) Rprintf("  while(gcm[ix] != lcm[iv]) :%s", (*debug >= 2) ? "\n" : " ");
				do { /* gcm[ix] != lcm[iv]  (after first loop) */
					/* long */ double dx;
					int gcmix = gcm[ix], lcmiv = lcm[iv];
					if (gcmix > lcmiv) {
						// If the next point of either the GCM or LCM is from the LCM, calculate the distance here.
						int gcmi1 = gcm[ix + 1];
						dx = (lcmiv - gcmi1 + 1) -
								(/*(long double)*/ x[lcmiv-1] - x[gcmi1-1]) * (gcmix - gcmi1)/(x[gcmix-1] - x[gcmi1-1]);
						++iv;
						if (dx >= d) {
							d = dx;
							ig = ix + 1;
							ih = iv - 1;
							//if(*debug >= 2) Rprintf(" L(%d,%d)", ig,ih);
						}
					}
					else {
						// If the next point of either the GCM or LCM is from the GCM, calculate the distance here.
						int lcmiv1 = lcm[iv - 1];
						/* Fix by Yong Lu {symmetric to above!}; original Fortran: only ")" misplaced! :*/
						dx = (/*(long double)*/x[gcmix-1] - x[lcmiv1-1]) * (lcmiv - lcmiv1) / (x[lcmiv-1] - x[lcmiv1-1])- (gcmix - lcmiv1 - 1);
						--ix;
						if (dx >= d) {
							d = dx;
							ig = ix + 1;
							ih = iv;
							//if(*debug >= 2) Rprintf(" G(%d,%d)", ig,ih);
						}
					}
					if (ix < 1)	ix = 1;
					if (iv > l_lcm)	iv = l_lcm;
					//if(*debug) {
					//	if(*debug >= 2) Rprintf(" --> ix = %d, iv = %d\n", ix,iv);
					//	else Rprintf(".");
					//}
				} while (gcm[ix] != lcm[iv]);
				//if(*debug && *debug < 2) Rprintf("\n");
			}
			else { /* l_gcm or l_lcm == 2 */
				d = 0; // d = (*min_is_0) ? 0. : 1.;
				//if(*debug) Rprintf("  ** (l_lcm,l_gcm) = (%d,%d) ==> d := %g\n", l_lcm, l_gcm, (double)d);
			}

			if (d < dip) break; // goto L_END;

			// Calculate the DIPs for the current LOW and HIGH.
			//if(*debug) Rprintf("  calculating dip ..");

			//int j_best, j_l = -1, j_u = -1;

			/* The DIP for the convex minorant. */
			dip_l = 0.;
			for (int j = ig; j < l_gcm; ++j) {
				double max_t = 1.;
				int /*j_ = -1,*/ jb = gcm[j + 1], je = gcm[j];
				if (je - jb > 1 && x[je-1] != x[jb-1]) {
					double C = (je - jb) / (x[je-1] - x[jb-1]);
					for (int jj = jb; jj <= je; ++jj) {
						double t = (jj - jb + 1) - (x[jj-1] - x[jb-1]) * C;
						if (max_t < t) {
							max_t = t; //j_ = jj;
						}
					}
				}
				if (dip_l < max_t) {
					dip_l = max_t; //j_l = j_;
				}
			}

			/* The DIP for the concave majorant. */
			dip_u = 0.;
			for (int j = ih; j < l_lcm; ++j) {
				double max_t = 1.;
				int /*j_ = -1,*/ jb = lcm[j], je = lcm[j + 1];
				if (je - jb > 1 && x[je-1] != x[jb-1]) {
					double C = (je - jb) / (x[je-1] - x[jb-1]);
					for (int jj = jb; jj <= je; ++jj) {
						double t = (x[jj-1] - x[jb-1]) * C - (jj - jb - 1);
						if (max_t < t) {
							max_t = t; //j_ = jj;
						}
					}
				}
				if (dip_u < max_t) {
					dip_u = max_t; //j_u = j_;
				}
			}

			//if(*debug) Rprintf(" (dip_l, dip_u) = (%g, %g)", dip_l, dip_u);

			/* Determine the current maximum. */
			if(dip_u > dip_l) {
				dipnew = dip_u; //j_best = j_u;
			} else {
				dipnew = dip_l; //j_best = j_l;
			}
			if (dip < dipnew) {
				dip = dipnew;
				//if(*debug) Rprintf(" -> new larger dip %g (j_best = %d)\n", dipnew, j_best);
			}
			//else if(*debug) Rprintf("\n");

			/*--- The following if-clause is NECESSARY  (may loop infinitely otherwise)!
		      --- Martin Maechler, Statistics, ETH Zurich, July 30 1994 ---------- */
			if (low == gcm[ig] && high == lcm[ih]) {
				//if(*debug) Rprintf("No improvement in  low = %ld  nor  high = %ld --> END\n", low, high);
				break;
			} else {
				low  = gcm[ig];
				high = lcm[ih];	// goto LOOP_Start; /* Recycle */
			}
		}
		/*---------------------------------------------------------------------------*/

		//L_END:
		/* do this in the caller :
		 *   *xl = x[low];  *xu = x[high];
		 * rather return the (low, high) indices -- automagically via lo_hi[]  */
		dip /= (2*n);
		int nn = qDiptabN.length - 1, max_n = qDiptabN[nn], prn = qDiptab[0].length;
		double zz[] = new double[prn];
		if (n > max_n) {
			double sqn0 = Math.sqrt(max_n);
			for (int j = 0; j < prn; j++)
				zz[j] = sqn0 * qDiptab[nn][j];
		} else {
			int i_n = nn;
			while(n < qDiptabN[i_n]) i_n--;
			int i2 = i_n + 1,
				n0 = qDiptabN[i_n],
				n1 = qDiptabN[i2];
			double f_n = (n - n0) * 1.0/(n1 - n0), sqn0 =Math.sqrt(n0), sqn1 = Math.sqrt(n1);
			for (int j = 0; j < prn; j++) {
				double y0 = sqn0 * qDiptab[i_n][j];
				zz[j] = y0 + f_n * (sqn1 * qDiptab[i2][j] - y0);
			}
		}
		double pval = 1 - linear(Math.sqrt(n) * dip, zz, qDiptabPr, 0, 1);
		return new double[] {dip, pval, low-1, high-1};
	}
	
	/**
	 * Linear approximation for diptest, also from jdistlib under gpl
	 * @param v
	 * @param x
	 * @param y
	 * @param lo
	 * @param hi
	 * @return Approximated value
	 */
	public static final double linear(double v, double[] x, double[] y, double lo, double hi) {
		int
			left = 0,
			right = x.length - 1;
		if (v < x[left])
			return lo;
		if (v > x[right])
			return hi;
		while(left < right - 1) {
			int mid = (left + right)/2;
			if(v < x[mid]) right = mid; else left = mid;
		}
		if(v == x[right])
			return y[right];
		if(v == x[left])
			return y[left];
		return v = y[left] + (y[right] - y[left]) * ((v - x[left])/(x[right] - x[left]));
	}
	
	private static final double[][] qDiptab = {
			{0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.125,0.132559548782689,0.157497369040235,0.187401878807559,0.20726978858736,0.223755804629222,0.231796258864192,0.237263743826779,0.241992892688593,0.244369839049632,0.245966625504691,0.247439597233262,0.248230659656638,0.248754269146416,0.249302039974259,0.249459652323225,0.24974836247845},
			{0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.108720593576329,0.121563798026414,0.134318918697053,0.147298798976252,0.161085025702604,0.176811998476076,0.186391796027944,0.19361253363045,0.196509139798845,0.198159967287576,0.199244279362433,0.199617527406166,0.199800941499028,0.199917081834271,0.199959029093075,0.199978395376082,0.199993151405815,0.199995525025673,0.199999835639211},
			{0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0833333333333333,0.0924514470941933,0.103913431059949,0.113885220640212,0.123071187137781,0.13186973390253,0.140564796497941,0.14941924112913,0.159137064572627,0.164769608513302,0.179176547392782,0.191862827995563,0.202101971042968,0.213015781111186,0.219518627282415,0.224339047394446,0.229449332154241,0.232714530449602,0.236548128358969,0.2390887911995,0.240103566436295,0.244672883617768},
			{0.0714285714285714,0.0714285714285714,0.0714285714285714,0.0725717816250742,0.0817315478539489,0.0940590181922527,0.103244490800871,0.110964599995697,0.117807846504335,0.124216086833531,0.130409013968317,0.136639642123068,0.144240669035124,0.159903395678336,0.175196553271223,0.184118659121501,0.191014396174306,0.198216795232182,0.202341010748261,0.205377566346832,0.208306562526874,0.209866047852379,0.210967576933451,0.212233348558702,0.212661038312506,0.21353618608817},
			{0.0625,0.0625,0.0656911994503283,0.0738651136071762,0.0820045917762512,0.0922700601131892,0.0996737189599363,0.105733531802737,0.111035129847705,0.115920055749988,0.120561479262465,0.125558759034845,0.141841067033899,0.153978303998561,0.16597856724751,0.172988528276759,0.179010413496374,0.186504388711178,0.19448404115794,0.200864297005026,0.208849997050229,0.212556040406219,0.217149174137299,0.221700076404503,0.225000835357532,0.233772919687683},
			{0.0555555555555556,0.0613018090298924,0.0658615858179315,0.0732651142535317,0.0803941629593475,0.0890432420913848,0.0950811420297928,0.0999380897811046,0.104153560075868,0.108007802361932,0.112512617124951,0.122915033480817,0.136412639387084,0.146603784954019,0.157084065653166,0.164164643657217,0.172821674582338,0.182555283567818,0.188658833121906,0.194089120768246,0.19915700809389,0.202881598436558,0.205979795735129,0.21054115498898,0.21180033095039,0.215379914317625},
			{0.05,0.0610132555623269,0.0651627333214016,0.0718321619656165,0.077966212182459,0.0852835359834564,0.0903204173707099,0.0943334983745117,0.0977817630384725,0.102180866696628,0.109960948142951,0.118844767211587,0.130462149644819,0.139611395137099,0.150961728882481,0.159684158858235,0.16719524735674,0.175419540856082,0.180611195797351,0.185286416050396,0.191203083905044,0.195805159339184,0.20029398089673,0.205651089646219,0.209682048785853,0.221530282182963},
			{0.0341378172277919,0.0546284208048975,0.0572191260231815,0.0610087367689692,0.0642657137330444,0.0692234107989591,0.0745462114365167,0.0792030878981762,0.083621033469191,0.0881198482202905,0.093124666680253,0.0996694393390689,0.110087496900906,0.118760769203664,0.128890475210055,0.13598356863636,0.142452483681277,0.150172816530742,0.155456133696328,0.160896499106958,0.166979407946248,0.17111793515551,0.175900505704432,0.181856676013166,0.185743454151004,0.192240563330562},
			{0.033718563622065,0.0474333740698401,0.0490891387627092,0.052719998201553,0.0567795509056742,0.0620134674468181,0.0660163872069048,0.0696506075066401,0.0733437740592714,0.0776460662880254,0.0824558407118372,0.088344627001737,0.0972346018122903,0.105130218270636,0.114309704281253,0.120624043335821,0.126552378036739,0.13360135382395,0.138569903791767,0.14336916123968,0.148940116394883,0.152832538183622,0.156010163618971,0.161319225839345,0.165568255916749,0.175834459522789},
			{0.0262674485075642,0.0395871890405749,0.0414574606741673,0.0444462614069956,0.0473998525042686,0.0516677370374349,0.0551037519001622,0.058265005347493,0.0614510857304343,0.0649164408053978,0.0689178762425442,0.0739249074078291,0.0814791379390127,0.0881689143126666,0.0960564383013644,0.101478558893837,0.10650487144103,0.112724636524262,0.117164140184417,0.121425859908987,0.126733051889401,0.131198578897542,0.133691739483444,0.137831637950694,0.141557509624351,0.163833046059817},
			{0.0218544781364545,0.0314400501999916,0.0329008160470834,0.0353023819040016,0.0377279973102482,0.0410699984399582,0.0437704598622665,0.0462925642671299,0.048851155289608,0.0516145897865757,0.0548121932066019,0.0588230482851366,0.0649136324046767,0.0702737877191269,0.0767095886079179,0.0811998415355918,0.0852854646662134,0.0904847827490294,0.0940930106666244,0.0974904344916743,0.102284204283997,0.104680624334611,0.107496694235039,0.11140887547015,0.113536607717411,0.117886716865312},
			{0.0164852597438403,0.022831985803043,0.0238917486442849,0.0256559537977579,0.0273987414570948,0.0298109370830153,0.0317771496530253,0.0336073821590387,0.0354621760592113,0.0374805844550272,0.0398046179116599,0.0427283846799166,0.047152783315718,0.0511279442868827,0.0558022052195208,0.059024132304226,0.0620425065165146,0.0658016011466099,0.0684479731118028,0.0709169443994193,0.0741183486081263,0.0762579402903838,0.0785735967934979,0.0813458356889133,0.0832963013755522,0.0926780423096737},
			{0.0111236388849688,0.0165017735429825,0.0172594157992489,0.0185259426032926,0.0197917612637521,0.0215233745778454,0.0229259769870428,0.024243848341112,0.025584358256487,0.0270252129816288,0.0286920262150517,0.0308006766341406,0.0339967814293504,0.0368418413878307,0.0402729850316397,0.0426864799777448,0.044958959158761,0.0477643873749449,0.0497198001867437,0.0516114611801451,0.0540543978864652,0.0558704526182638,0.0573877056330228,0.0593365901653878,0.0607646310473911,0.0705309107882395},
			{0.00755488597576196,0.0106403461127515,0.0111255573208294,0.0119353655328931,0.0127411306411808,0.0138524542751814,0.0147536004288476,0.0155963185751048,0.0164519238025286,0.017383057902553,0.0184503949887735,0.0198162679782071,0.0218781313182203,0.0237294742633411,0.025919578977657,0.0274518022761997,0.0288986369564301,0.0306813505050163,0.0320170996823189,0.0332452747332959,0.0348335698576168,0.0359832389317461,0.0369051995840645,0.0387221159256424,0.03993025905765,0.0431448163617178},
			{0.00541658127872122,0.00760286745300187,0.00794987834644799,0.0085216518343594,0.00909775605533253,0.00988924521014078,0.0105309297090482,0.0111322726797384,0.0117439009052552,0.012405033293814,0.0131684179320803,0.0141377942603047,0.0156148055023058,0.0169343970067564,0.018513067368104,0.0196080260483234,0.0206489568587364,0.0219285176765082,0.0228689168972669,0.023738710122235,0.0248334158891432,0.0256126573433596,0.0265491336936829,0.027578430100536,0.0284430733108,0.0313640941982108},
			{0.00390439997450557,0.00541664181796583,0.00566171386252323,0.00607120971135229,0.0064762535755248,0.00703573098590029,0.00749421254589299,0.00792087889601733,0.00835573724768006,0.00882439333812351,0.00936785820717061,0.01005604603884,0.0111019116837591,0.0120380990328341,0.0131721010552576,0.0139655122281969,0.0146889122204488,0.0156076779647454,0.0162685615996248,0.0168874937789415,0.0176505093388153,0.0181944265400504,0.0186226037818523,0.0193001796565433,0.0196241518040617,0.0213081254074584},
			{0.00245657785440433,0.00344809282233326,0.00360473943713036,0.00386326548010849,0.00412089506752692,0.00447640050137479,0.00476555693102276,0.00503704029750072,0.00531239247408213,0.00560929919359959,0.00595352728377949,0.00639092280563517,0.00705566126234625,0.0076506368153935,0.00836821687047215,0.00886357892854914,0.00934162787186159,0.00993218636324029,0.0103498795291629,0.0107780907076862,0.0113184316868283,0.0117329446468571,0.0119995948968375,0.0124410052027886,0.0129467396733128,0.014396063834027},
			{0.00174954269199566,0.00244595133885302,0.00255710802275612,0.00273990955227265,0.0029225480567908,0.00317374638422465,0.00338072258533527,0.00357243876535982,0.00376734715752209,0.00397885007249132,0.00422430013176233,0.00453437508148542,0.00500178808402368,0.00542372242836395,0.00592656681022859,0.00628034732880374,0.00661030641550873,0.00702254699967648,0.00731822628156458,0.0076065423418208,0.00795640367207482,0.0082270524584354,0.00852240989786251,0.00892863905540303,0.00913853933000213,0.00952234579566773},
			{0.00119458814106091,0.00173435346896287,0.00181194434584681,0.00194259470485893,0.00207173719623868,0.00224993202086955,0.00239520831473419,0.00253036792824665,0.00266863168718114,0.0028181999035216,0.00299137548142077,0.00321024899920135,0.00354362220314155,0.00384330190244679,0.00420258799378253,0.00445774902155711,0.00469461513212743,0.00499416069129168,0.00520917757743218,0.00540396235924372,0.00564540201704594,0.00580460792299214,0.00599774739593151,0.00633099254378114,0.00656987109386762,0.00685829448046227},
			{0.000852415648011777,0.00122883479310665,0.00128469304457018,0.00137617650525553,0.00146751502006323,0.00159376453672466,0.00169668445506151,0.00179253418337906,0.00189061261635977,0.00199645471886179,0.00211929748381704,0.00227457698703581,0.00250999080890397,0.00272375073486223,0.00298072958568387,0.00315942194040388,0.0033273652798148,0.00353988965698579,0.00369400045486625,0.00383345715372182,0.00400793469634696,0.00414892737222885,0.0042839159079761,0.00441870104432879,0.00450818604569179,0.00513477467565583},
			{0.000644400053256997,0.000916872204484283,0.000957932946765532,0.00102641863872347,0.00109495154218002,0.00118904090369415,0.00126575197699874,0.00133750966361506,0.00141049709228472,0.00148936709298802,0.00158027541945626,0.00169651643860074,0.00187306184725826,0.00203178401610555,0.00222356097506054,0.00235782814777627,0.00248343580127067,0.00264210826339498,0.0027524322157581,0.0028608570740143,0.00298695044508003,0.00309340092038059,0.00319932767198801,0.00332688234611187,0.00339316094477355,0.00376331697005859}};
		private static final double[] qDiptabPr = {0,0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98,0.99,0.995,0.998,0.999,0.9995,0.9998,0.9999,0.99995,0.99998,0.99999,1};
		private static final int[] qDiptabN = {4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 40000, 72000};

}
