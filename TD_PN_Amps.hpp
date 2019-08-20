/*
 * TD_PN_Amps.hpp
 *
 *  Created on: Feb 21, 2019
 *      Author: blakemoore
 */

#ifndef TD_PN_AMPS_HPP_
#define TD_PN_AMPS_HPP_
#include <math.h>

static double Bpow(double a, double b){
	if(a > 1 + 1E-5 || a < 1 - 1E-5 ){
		return pow(a, b);
	} else {
		double ord0 = 1.;
		double ord1 = b*(a-1.);
		double ord2 = 1./2.*(b*b-b)*(a-1.)*(a-1.);
		double ord3 = 1./6.*(b - 2)*(b - 1)*b*(a-1)*(a-1)*(a-1);

		return ord0 + ord1 + ord2 + ord3;
	}
}

static double  P_0_C_2_C_2(double e, double u, double iota){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = 1./Bpow(1-ecosu, 2)*((1+C*C)*(ecosu*ecosu-ecosu-2*e*e+2));
	return term;
}
static double  P_0_C_2_S_2(double e, double u, double iota){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = 2./Bpow(1-ecosu, 2)*((1+C*C)*Bpow(1-e*e,1./2.)*esinu);
	return term;
}
static double  P_0(double e, double u, double iota){
	double S = sin(iota);
	double ecosu = e*cos(u);
	double term;

	term = S*S/(1-ecosu)*ecosu;
	return term;
}
static double  P_05_C_1_C_1(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double M = m1 + m2;
	double term;

	term = S/2*(m1-m2)/M*1./(1-ecosu)*((1-3*C*C)*esinu);
	return term;
}
static double  P_05_C_1_S_1(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = S/4*(m1-m2)/M*1./Bpow(1-ecosu,2)*(Bpow(1-e*e,1./2.)*(2*ecosu+5-(6*ecosu-1)*C*C));
	return term;
}
static double  P_05_C_3_C_3(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double S = sin(iota);
	double C = cos(iota);
	double term;
	double M = m1 + m2;

	term = -S/2*(m1-m2)/M*1./Bpow(1-ecosu,3)*((1+C*C)*esinu*(ecosu*ecosu-2*ecosu-4*e*e+5));
	return term;
}
static double  P_05_C_3_S_3(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = S/4*(m1-m2)/M*1./Bpow(1-ecosu,3)*((1+C*C)*Bpow(1-e*e,1./2.)*(6*ecosu*ecosu-7*ecosu-8*e*e+9));
	return term;
}
static double  P_1_C_2_C_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = 1./6.*1./Bpow(1-ecosu,3)*((-(9-5*eta)*pow(ecosu,3)+(18-10*eta)*ecosu*ecosu+((18-10*eta)*e*e+20-12*eta)*ecosu-(33+11*eta)*e*e-14+38*eta)
		+(-(3+13*eta)*pow(ecosu,3)+(6+26*eta)*ecosu*ecosu+((6+26*eta)*e*e+33-51*eta)*ecosu-(48-34*eta)*e*e+6-22*eta)*C*C+(1-3*eta)*(-6*pow(ecosu,3)
		+12*ecosu*ecosu+(12*e*e-13)*ecosu-9*e*e+4)*pow(C,4));
	return term;
}
static double  P_1_C_2_S_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = 1./6.*esinu/Bpow(1-ecosu,3)*(1./Bpow(1-e*e,1./2.)*(((18-10*eta)*e*e-6-2*eta)*ecosu-(39-7*eta)*e*e+27+5*eta)+1./Bpow(1-e*e,1./2.)*(((6+26*eta)*e*e
			+6-38*eta)*ecosu-(30+20*eta)*e*e+18+32*eta)*C*C+3*Bpow(1-e*e,1./2.)*(1-3*eta)*(-4*ecosu+5)*pow(C,4));
	return term;
}
static double  P_1_C_4_C_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/24*1./Bpow(1-ecosu,4)*((1-3*eta)*(1+C*C)*(-6*pow(ecosu,4)+18*pow(ecosu,3)+(48*e*e-61)*ecosu*ecosu+(65-69*e*e)*ecosu-48*pow(e,4)+117*e*e-64));
	return term;
}
static double  P_1_C_4_S_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/4*1./Bpow(1-ecosu,4)*((1+C*C)*(1-3*eta)*Bpow(1-e*e,1./2.)*esinu*(-4*ecosu*ecosu+9*ecosu+8*e*e-13));
	return term;
}
static double  P_1(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/24*1./Bpow(1-ecosu,3)*(((30-2*eta)*pow(ecosu,3)-(60-4*eta)*ecosu*ecosu-(57+5*eta)*ecosu+(87+3*eta)*e*e)+
			((18-54*eta)*pow(ecosu,3)-(36-108*eta)*ecosu*ecosu+(3-9*eta)*ecosu+(15-45*eta)*e*e)*C*C);
	return term;
}
static double  P_15_C_1_C_1(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = -S/48*(m1-m2)/M*esinu/Bpow(1-ecosu,4)*((48*pow(ecosu,3)-144*ecosu*ecosu+(33+22*eta)*ecosu-(336+12*eta)*e*e+399-10*eta)
			+(-(108+72*eta)*pow(ecosu,3)+(324+216*eta)*ecosu*ecosu-(12+240*eta)*ecosu-(144-36*eta)*e*e-60+60*eta)*C*C
			+5*(1-ecosu)*(1-2*eta)*(12*ecosu*ecosu-24*ecosu+5)*C*C*C*C);
	return term;
}
static double  P_15_C_1_S_1(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = S/96*(m1-m2)/M*Bpow(1-ecosu,-4)*(Bpow(1-e*e,-1./2.)*((48+48*eta-96*e*e)*pow(ecosu,3)-(24+48*eta-72*e*e*eta)*ecosu*ecosu
			+((396+104*eta)*e*e-204-296*eta)*ecosu+(123-198*eta)*pow(e,4)-(546-220*eta)*e*e+303+98*eta)+Bpow(1-e*e,-1./2.)*(((216+144*eta)*e*e-72-288*eta)*pow(ecosu,3)
					-((756+240*eta)*e*e-444-552*eta)*ecosu*ecosu-((144+96*eta)*e*e-336+96*eta)*ecosu+(720+144*eta)*pow(e,4)-(756+96*eta)*e*e+12-24*eta)*C*C
					-Bpow(1-e*e,1./2.)*(1-2*eta)*(120*pow(ecosu,3)-276*ecosu*ecosu+52*ecosu+105*e*e-1)*pow(C,4));
	return term;
}
static double  P_15_C_3_C_3(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = S/96*(m1-m2)/M*esinu*Bpow(1-ecosu,-4)*(((108-24*eta)*pow(ecosu,3)-(324-72*eta)*ecosu*ecosu-((432-96*eta)*e*e-143-274*eta)*ecosu+(1176-72*eta)*e*e
			-671-346*eta)+((72+48*eta)*pow(ecosu,3)-(216+144*eta)*ecosu*ecosu-((288+192*eta)*e*e+88-736*eta)*ecosu+(1152-24*eta)*e*e-632-424*eta)*C*C
			+5*(1-2*eta)*(12*pow(ecosu,3)-36*ecosu*ecosu+(77-48*e*e)*ecosu+72*e*e-77)*pow(C,4));
	return term;
}
static double  P_15_C_3_S_3(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = S/64*(m1-m2)/M*Bpow(1-ecosu,-4)*(Bpow(1-e*e,-1./2.)*(((216-48*eta)*e*e-120-48*eta)*pow(ecosu,3)-((644-120*eta)*e*e-436-88*eta)*ecosu*ecosu
			-((288-64*eta)*pow(e,4)-(76+232*eta)*e*e-340+424*eta)*ecosu+(719+2*eta)*pow(e,4)-(510+436*eta)*e*e-225+450*eta)
			+Bpow(1-e*e,-1/2.)*(((144+96*eta)*e*e-48-192*eta)*pow(ecosu,3)-((504+160*eta)*e*e-296-368*eta)*ecosu*ecosu-((192+128*eta)*pow(e,4)
					+(96-576*eta)*e*e-416+576*eta)*ecosu+(928-416*eta)*pow(e,4)-(1016-576*eta)*e*e+72-144*eta)*C*C-Bpow(1-e*e,1./2.)*(1-2*eta)*(120*pow(ecosu,3)
	-276*ecosu*ecosu-(160*e*e-212)*ecosu+185*e*e-81)*pow(C,4));
	return term;
}
static double  P_15_C_5_C_5(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = pow(S,3)/96*(m1-m2)/M*esinu*Bpow(1-ecosu,-5)*(1+C*C)*(1-2*eta)*(12*pow(ecosu,4)-48*pow(ecosu,3)-(144*e*e-209)*ecosu*ecosu
			+(360*e*e-394)*ecosu+192*pow(e,4)-600*e*e+413);
	return term;
}
static double  P_15_C_5_S_5(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = -pow(S,3)/192*(m1-m2)/M*Bpow(1-e*e,1./2.)*Bpow(1-ecosu,-5)*(1+C*C)*(1-2*eta)*(120*pow(ecosu,4)-396*pow(ecosu,3)
	-(480*e*e-808)*ecosu*ecosu+(825*e*e-773)*ecosu+384*pow(e,4)-1113*e*e+625);
	return term;
}
static double  P_2_C_2_C_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = 1./5760.*Bpow(1-ecosu,-5)*(2880*Bpow(1-e*e,-1./2.)*(1+C*C)*Bpow(1-ecosu,2)*(5-2*eta)*(-2*pow(ecosu,3)+7*ecosu*ecosu+(4*e*e+3)*ecosu
			-16*e*e+4)+1/(1-e*e)*(-(1-e*e)*(7560-15720*eta-600*eta*eta)*pow(ecosu,5)+(1-e*e)*(30240-62880*eta-2400*eta*eta)*pow(ecosu,4)
					+(-(15120-31440*eta-1200*eta*eta)*pow(e,4)+
			+(7518-66070*eta+5670*eta*eta)*e*e-65838+82150*eta-6870*eta*eta)*pow(ecosu,3)+((28980+15420*eta-35580*eta*eta)*pow(e,4)
			+(64686-152550*eta-6090*eta*eta)*e*e+126654-5430*eta+41670*eta*eta)*ecosu*ecosu+((99585+20955*eta-19755*eta*eta)*pow(e,4)
			-(293136-38160*eta-163200*eta*eta)*e*e-26769+83445*eta-143445*eta*eta)*ecosu+(10665-115245*eta+63405*eta*eta)*pow(e,6)
			-(145440-277920*eta+136080*eta*eta)*pow(e,4)+(275607-212435*eta+25635*eta*eta)*e*e-67392+2240*eta+47040*eta*eta)
			+1/(1-e*e)*(((1-e*e)*(6120-21960*eta-11640*eta*eta))*pow(ecosu,5)-(1-e*e)*(24480-87840*eta-46560*eta*eta)*pow(ecosu,4)
			+((12240-43920*eta-23280*eta*eta)*pow(e,4)-(75006-129750*eta-181050*eta*eta)*e*e-10674-38310*eta-157770*eta*eta)*pow(ecosu,3)
			+(-(16020-113700*eta-56220*eta*eta)*pow(e,4)+(108018-202650*eta-213270*eta*eta)*e*e+128322-53610*eta+157050*eta*eta)*ecosu*ecosu
			+((67455+118965*eta+7995*eta*eta)*pow(e,4)-(202608+171600*eta+89760*eta*eta)*e*e-85167+195195*eta+81765*eta*eta)*ecosu
			+(12375-88515*eta-64125*eta*eta)*pow(e,6)-(100800-76800*eta-151440*eta*eta)*pow(e,4)+(188361+44835*eta-35475*eta*eta)*e*e
			-26496-80640*eta-51840*eta*eta)*C*C+5*(-(1800-2856*eta-7128*eta*eta)*pow(ecosu,5)+(7200-11424*eta-28512*eta*eta)*pow(ecosu,4)
			+((3600-5712*eta-14256*eta*eta)*e*e-4746-7342*eta+61614*eta*eta)*pow(ecosu,3)+(-(18684-43500*eta-34452*eta*eta)*e*e-6
			+23358*eta-67518*eta*eta)*ecosu*ecosu+((12909-32529*eta-16815*eta*eta)*e*e+8109-30609*eta+19281*eta*eta)*ecosu+(9045
			-24345*eta-7335*eta*eta)*pow(e,4)-(15915-43431*eta-11289*eta*eta)*e*e+288-1184*eta+672*eta*eta)*pow(C,4)
			+15*(1-5*eta+5*eta*eta)*(-360*pow(ecosu,5)+1440*pow(ecosu,4)+(720*e*e-2418)*pow(ecosu,3)-(2124*e*e-2178)*ecosu*ecosu
			+(1553*e*e-527)*ecosu+345*pow(e,4)-839*e*e+32)*pow(C,6));
			return term;
}

// checked to here

static double  P_2_C_2_S_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = 1./192*esinu*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-1./2.)*((-(504-1048*eta-40*eta*eta)*pow(e,4)+(144-592*eta+48*eta*eta)*e*e-18080+216*eta
			+104*eta*eta)*pow(ecosu,3)+((1584-200*eta-448*eta*eta)*pow(e,4)-(1536-1568*eta-416*eta*eta)*e*e+4272-784*eta-544*eta*eta)*ecosu*ecosu
			+((1956+908*eta-92*eta*eta)*pow(e,4)-(4584-784*eta-544*eta*eta)*ecosu*ecosu+((1956+908*eta-92*eta*eta)*pow(e,4)-(4584-1736*eta-760*eta*eta)*e*e
			-1692-628*eta-92*eta*eta)*ecosu+(1359-6699*eta+171*eta*eta)*pow(e,6)-(7113-20941*eta+13*eta*eta)*pow(e,4)+(10053-22809*eta-711*eta*eta)*pow(e,4)
			-(1296-2896*eta-2832*eta*eta)*e*e-552-760*eta-1864*eta*eta)*pow(ecosu,3)+(-(1536-5824*eta-2000*eta*eta)*pow(e,4)+(3648-11456*eta-7648*eta*eta)*e*e
			+2208+3616*eta+5072*eta*eta)*ecosu*ecosu+((1692+2740*eta-3012*eta*eta)*pow(e,4)-(3096+5768*eta-9480*eta*eta)*e*e-2916+5044*eta-5892*eta*eta)*ecosu
			-(1167-987*eta-1605*eta*eta)*pow(e,6)+(2937-10061*eta-3027*eta*eta)*pow(e,4)-(2757-17289*eta-151*eta*eta)*e*e+2427-8887*eta+1079*eta*eta)*C*C
			+Bpow(1-e*e,-1./2.)*(((600-952*eta-2376*eta*eta)*e*e-216-584*eta+3528*eta*eta)*pow(ecosu,3)+(-(2688-5440*eta-7296*eta*eta)*e*e+1440-448*eta
			-11040*eta*eta)*ecosu*ecosu+((1684-2596*eta-6940*eta*eta)*e*e-340-2780*eta+10972*eta*eta)*ecosu+(975-3195*eta+795*eta*eta)*pow(e,4)-(1546-4498*eta
			-430*eta*eta)*e*e+91+617*eta-2665*eta*eta)*pow(C,4)+5*(1-5*eta+5*eta*eta)*Bpow(1-e*e,1./2.)*(-72*pow(ecosu,3)+240*ecosu*ecosu-220*ecosu+3*e*e+49)*pow(C,6))
			+(1+C*C)*(5-2*eta)*Bpow(1-e*e,-1)*esinu*Bpow(1-ecosu,4)*((2*e*e+1)*ecosu-8*e*e+5);
	return term;
}
static double  P_2_C_4_C_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/2880*Bpow(1-ecosu,-5)*(((2160-6960*eta+720*eta*eta)*pow(ecosu,5)-(8640-27840*eta+2880*eta*eta)*pow(ecosu,4)+(-(17280-55680*eta+5760*eta*eta)*e*e+14238-37370*eta
			-25350*eta*eta)*pow(ecosu,3)+((63000-199800*eta+20520*eta*eta)*e*e-31314+78150*eta+59850*eta*eta)*ecosu*ecosu+((17280-55680*eta+5760*eta*eta)*pow(e,4)
			-(24105-60435*eta-49125*eta*eta)*e*e-23661+107055*eta-114375*eta*eta)*ecosu-(49905-148875*eta+7155*eta*eta)*pow(e,4)+(43635-102705*eta-61095*eta*eta)*e*e+14592
			-75520*eta+80640*eta*eta)+((1800-5160*eta-1080*eta*eta)*pow(ecosu,5)-(7200-20640*eta-4320*eta*eta)*pow(ecosu,4)+(-(14400-41280*eta-8640*eta*eta)*e*e+9660-14480*eta
			-48240*eta*eta)*pow(ecosu,3)+((58140-175500*eta-3780*eta*eta)*e*e-26400+53580*eta+84420*eta*eta)*ecosu*ecosu+((14400-41280*eta-8640*eta*eta)*pow(e,4)-(21780
			-48810*eta-60750*eta*eta)*e*e-22080+99150*eta-106470*eta*eta)*ecosu-(62460-211650*eta+69930*eta*eta)*pow(e,4)+(74160-255330*eta+91530*eta*eta)*e*e-3840+16640*eta
			-11520*eta*eta)*C*C+3*(1-5*eta+5*eta*eta)*(360*pow(ecosu,5)-1440*pow(ecosu,4)-(2880*e*e-4578*pow(ecosu,3)+(7740*e*e-7794)*ecosu*ecosu+(2880*pow(e,4)
			-8885*e*e+4979)*ecosu-4245*pow(e,4)+6755*e*e-2048)*pow(C,4)));
	return term;
}
static double  P_2_C_4_S_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/48*esinu*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-1./2.)*((-(144-464*eta+48*eta*eta)*e*e+96-272*eta-96*eta*eta)*pow(ecosu,3)+((558-1790*eta+198*eta*eta)*e*e
			-402+1166*eta+270*eta*eta)*ecosu*ecosu+((288-928*eta+96*eta*eta)*pow(e,4)-(680-1936*eta-632*eta*eta)*e*e+224-336*eta-1232*eta*eta)*ecosu-(879-2665*eta
			+93*eta*eta)*pow(e,4)+(1448-4084*eta-788*eta*eta)*e*e-509+1179*eta+1061*eta*eta)+Bpow(1-e*e,-1./2.)*((-(120-344*eta-72*eta*eta)*e*e+72-152*eta-216*eta*eta)*pow(ecosu,3)
			+((486-1430*eta-162*eta*eta)*e*e-330+806*eta+630*eta*eta)*ecosu*ecosu+((240-688*eta-144*eta*eta)*pow(e,4)-(540-1236*eta-1332*eta*eta)*e*e+132+124*eta
			-1692*eta*eta)*ecosu-(936-2950*eta+378*eta*eta)*pow(e,4)+(1566-4674*eta-198*eta*eta)*e*e-570+1484*eta+756*eta*eta)*C*C+Bpow(1-e*e,1./2.)*(1-5*eta
			+5*eta*eta)*(72*pow(ecosu,3)-240*ecosu*ecosu-(144*e*e-364)*ecosu+249*e*e-301)*pow(C,4));
	return term;
}
static double  P_2_C_6_C_6(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = pow(S,4)/1920.*(1+C*C)*(1-5*eta+5*eta*eta)*Bpow(1-ecosu,-6)*(120*pow(ecosu,6)-600*pow(ecosu,5)-(2160*e*e-3206)*pow(ecosu,4)+(7860*e*e-8444)*pow(ecosu,3)
			+(5760*pow(e,4)-19135*e*e+13051)*ecosu*ecosu-(11475*pow(e,4)-23240*e*e+11269)*ecosu-3840*pow(e,6)+17235*pow(e,4)-21325*e*e+7776);
	return term;
}
static double  P_2_C_6_S_6(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = Bpow(1-e*e,1./2.)/192*pow(S,4)*(1+C*C)*(1-5*eta+5*eta*eta)*esinu*Bpow(1-ecosu,-6)*(72*pow(ecosu,4)-312*pow(ecosu,3)-(384*e*e-844)*ecosu*ecosu
			+(1053*e*e-1325)*ecosu+384*pow(e,4)-1437*e*e+1105);
	return term;
}
static double  P_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = S*S/192*Bpow(1-ecosu,-5)*(((120-120*eta-8*eta*eta)*pow(ecosu,5)+(3360-2608*eta+416*eta*eta)*pow(ecosu,4)-(11058-7718*eta+1086*eta*eta)*pow(ecosu,3)
			+((168+1240*eta+152*eta*eta)*e*e+13710-8746*eta+738*eta*eta)*ecosu*ecosu+(-(3441-139*eta-277*eta*eta)*e*e-5181+2751*eta-423*eta*eta)*ecosu
			+(135-477*eta-363*eta*eta)*pow(e,4)+(3003-425*eta+297*eta*eta)*e*e-816+528*eta)+2*((180-516*eta-108*eta*eta)*pow(ecosu,5)
			-(720-2064*eta-432*eta*eta)*pow(ecosu,4)+(678-1928*eta-504*eta*eta)*pow(ecosu,3)+((150-430*eta-90*eta*eta)*e*e+336-1010*eta+90*eta*eta)*ecosu*ecosu
			+(-(930-2705*eta-315*eta*eta)*e*e-96+283*eta+9*eta*eta)*ecosu+(378-1107*eta-81*eta*eta)*pow(e,4)+(24-61*eta-63*eta*eta)*e*e)*C*C+(1-5*eta
			+5*eta*eta)*(120*pow(ecosu,5)-480*pow(ecosu,4)+566*pow(ecosu,3)+(84*e*e-102)*ecosu*ecosu-(343*e*e-1)*ecosu+105*pow(e,4)+49*e*e)*pow(C,4))
			+S*S/12*Bpow(1-ecosu,-2)*(Bpow(1-e*e,-1)*(((1-e*e)*(240-193*eta+24*eta*eta))*ecosu+51-33*eta)+ecosu*Bpow(1-e*e,-1./2.)*(-(60-24*eta)*ecosu+150-60*eta));
	return term;
}
static double  X_0_C_2_C_2(double e, double u, double iota){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = -4*C*Bpow(1-e*e,1./2.)*Bpow(1-ecosu,-2)*esinu;
	return term;
}
static double  X_0_C_2_S_2(double e, double u, double iota){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = 2*C*Bpow(1-ecosu,-2)*(ecosu*ecosu-ecosu-2*(e*e-1));
	return term;
}
static double  X_05_C_1_C_1(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*S/2*(m1-m2)/M*Bpow(1-e*e,1./2.)*Bpow(1-ecosu,-2)*(2*ecosu-3);
	return term;
}
static double  X_05_C_1_S_1(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = -C*S*(m1-m2)/M*Bpow(1-ecosu,-1)*esinu;
	return term;
}
static double  X_05_C_3_C_3(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = -C*S/2.*(m1-m2)/M*Bpow(1-ecosu,-3)*(Bpow(1-e*e,1./2.)*(6*ecosu*ecosu-7*ecosu-8*e*e+9));
	return term;
}
static double  X_05_C_3_S_3(double e, double u, double iota, double m1, double m2){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = -C*S*(m1-m2)/M*Bpow(1-ecosu,-3)*(esinu*(ecosu*ecosu-2*ecosu-4*e*e+5));
	return term;
}
static double  X_1_C_2_C_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = C/3.*esinu*Bpow(1-ecosu,-3)*(Bpow(1-e*e,-1./2.)*((-(12+8*eta)*e*e+20*eta)*ecosu+(33+11*eta)*e*e-21-23*eta)+3*Bpow(1-e*e,1./2.)*(1-3*eta)*(2*ecosu-3)*C*C);
	return term;
}
static double  X_1_C_2_S_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = C/6*Bpow(1-ecosu,-3)*((-(12+8*eta)*pow(ecosu,3)+(24+16*eta)*ecosu*ecosu+((24+16*eta)*e*e+53-63*eta)*ecosu-(69+13*eta)*e*e-20+52*eta)+(1-3*eta)*(-6*pow(ecosu,3)
	+12*ecosu*ecosu-(13-12*e*e)*ecosu-21*e*e+16)*C*C);
	return term;
}
static double  X_1_C_4_C_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S*S/2.*esinu*Bpow(1-ecosu,-4)*Bpow(1-e*e,1./2.)*(1-3*eta)*(4*ecosu*ecosu-9*ecosu-8*e*e+13);
	return term;
}
static double  X_1_C_4_S_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S/12*Bpow(1-ecosu,4)*(1-3*eta)*(-6*pow(ecosu,4)+18*pow(ecosu,3)+(48*e*e-61)*ecosu*ecosu-(69*e*e-65)*ecosu-48*pow(e,4)*117*e*e-64);
	return term;
}
static double  X_1(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S*S/2*esinu*Bpow(1-ecosu,-3)*Bpow(1-e*e,1./2.)*(1-3*eta);
	return term;
}
static double  X_15_C_1_C_1(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*S/48*(m1-m2)/M*Bpow(1-ecosu,-4)*(Bpow(1-e*e,-1./2.)*(-(96*e*e-48-48*eta)*pow(ecosu,3)+((432-24*eta)*e*e-264-144*eta)*ecosu*ecosu
			-((84+88*eta)*e*e+108-280*eta)*ecosu-(477-138*eta)*pow(e,4)+(702-164*eta)*e*e-153-46*eta)+Bpow(1-e*e,1./2.)*(1-2*eta)*(24*pow(ecosu,3)
			-84*ecosu*ecosu+68*ecosu-3*e*e-5)*C*C);
	return term;
}
static double  X_15_C_1_S_1(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*S/24*(m1-m2)/M*esinu*Bpow(1-ecosu,-4)*((48*pow(ecosu,3)-144*ecosu*ecosu+(33+22*eta)*ecosu+216+36*eta*e*e-153-58*eta)
			+(1-2*eta)*(12*pow(ecosu,3)-36*ecosu*ecosu+29*ecosu+24*e*e-29)*C*C);
	return term;
}
static double  X_15_C_3_C_3(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*S/32*(m1-m2)/M*Bpow(1-ecosu,-4)*(Bpow(1-e*e,-1./2.)*(((168+48*eta)*e*e-72-144*eta)*pow(ecosu,3)+(-(540+88*eta)*e*e+332+296*eta)*ecosu*ecosu
			+(-(224+64*eta)*pow(e,4)-(60-504*eta)*e*e+412-568*eta)*ecosu+(725-10*eta)*pow(e,4)-(570+316*eta)*e*e-171+342*eta)+Bpow(1-e*e,1./2.)*(1-2*eta)*(-72*pow(ecosu,3)
			+172*ecosu*ecosu+(96*e*e-140)*ecosu-191*e*e+135)*C*C);
	return term;
}
static double  X_15_C_3_S_3(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*S/48*(m1-m2)/M*esinu*Bpow(1-ecosu,4)*(((84+24*eta)*pow(ecosu,3)-(253+72*eta)*ecosu*ecosu-((336+96*eta)*e*e+11-582*eta)*ecosu+(1080+120*eta)*e*e-565-558*eta)
			+3*(1-2*eta)*(12*pow(ecosu,3)-36*ecosu*ecosu-(48*e*e-77)*ecosu+88*e*e-93)*C*C);
	return term;
}
static double  X_15_C_5_C_5(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*pow(S,3)/96*(m1-m2)/M*Bpow(1-ecosu,5)*(1-2*eta)*(Bpow(1-e*e,1./2.)*(120*pow(ecosu,4)-396*pow(ecosu,3)-(480*e*e-808)*ecosu*ecosu
			+(825*e*e-773)*ecosu+384*pow(e,4)-1113*e*e+625));
	return term;
}
static double  X_15_C_5_S_5(double e, double u, double iota, double m1, double m2, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;
	double M = m1 + m2;

	term = C*pow(S,3)/48*(m1-m2)/M*esinu*Bpow(1-ecosu,-5)*(1-2*eta)*(12*pow(ecosu,4)-48*pow(ecosu,3)-(144*e*e-209)*ecosu*ecosu+(360*e*e-394)*ecosu+192*pow(e,4)-600*e*e+413);
	return term;
}
static double  X_2_C_2_C_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double term;

	term = 2*C*esinu*Bpow(1-ecosu,-3)*((5-2*eta)/(1-e*e)*(-(2*e*e+1)*ecosu+8*e*e-5))+C/96*esinu*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-3./2.)*((-(24-568*eta-8*eta*eta)*pow(e,4)+(720-1872*eta
			-720*eta*eta)*e*e+744+632*eta+520*eta*eta)*pow(ecosu,3)+((336-3056*eta+352*eta*eta)*pow(e,4)-(1728-7840*eta-1504*eta*eta)*e*e-2928-2768*eta-1280*eta*eta)*ecosu*ecosu
			+(-(2652-988*eta-852*eta*eta)*pow(e,4)+(5400-3224*eta-4008*eta*eta)*e*e+1572+220*eta+2580*eta*eta)*ecosu-(885-5145*eta+297*eta*eta)*pow(e,6)+(4995-13935*eta
			-321*eta*eta)*pow(e,4)-(7047-12691*eta-2333*eta*eta)*e*e+1497-3229*eta-1523*eta*eta)+2./Bpow(1-e*e,1./2.)*((-(216-568*eta-264*eta*eta)*e*e-600+1208*eta
			+1848*eta*eta)*ecosu*ecosu+(-(868-2220*eta-1220*eta*eta)*e*e+484-684*eta-2372*eta*eta)*ecosu-(597-1737*eta-303*eta*eta)*pow(e,4)+(1342-3710*eta-1250*eta*eta)*e*e
			-601+1397*eta+1379*eta*eta)*C*C+Bpow(1-e*e,1./2.)*(1-5*eta+5*eta*eta)*(120*pow(ecosu,3)-432*ecosu*ecosu+484*ecosu+75*e*e-247)*pow(C,4));
	return term;
}
static double  X_2_C_2_S_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double term;

	term = C/2880*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-1./2.)*(((1-e*e)*(360-8520*eta-120*eta*eta))*pow(ecosu,5)-((1-e*e)*(1440-34080*eta-480*eta*eta))*pow(ecosu,4)
			+((720-17040*eta-240*eta*eta)*pow(e,4)-(43158-78910*eta-46290*eta*eta)*e*e-31002-14350*eta-46050*eta*eta)*pow(ecosu,3)+(-(16020-146340*eta
			+21540*eta*eta)*pow(e,4)+(119994-295890*eta-79710*eta*eta)*e*e+116346+6990*eta+101250*eta*eta)*ecosu*ecosu+((110835-10935*eta-18105*eta*eta)*pow(e,4)
			-(276384+12480*eta-145680*eta*eta)*e*e-54771+165975*eta-127575*eta*eta)*ecosu+(26955-161775*eta+45135*eta*eta)*pow(e,6)-(176400-366960*eta+95520*eta*eta)*pow(e,4)
			+(279333-230305*eta+23505*eta*eta)*e*e-56448-22400*eta+26880*eta*eta)+Bpow(1-e*e,-3/2)*((-(1-e*e)*(28800-11520*eta))*pow(ecosu,5)+((1-e*e)*(158400
			-63360*eta))*pow(ecosu,4)-((57600-23040)*eta*pow(e,4)-(244800-97920*eta)*e*e+187200-74880*eta)*pow(ecosu,3)+((345600-138240*eta)*pow(e,4)-(417600-167040*eta)*e*e
			+72000-28800*eta)*ecosu*ecosu+(-(518400-207360*eta)*pow(e,4)+(590400-236160*eta)*e*e-72000+28800*eta)*ecosu+(230400-92160*eta)*pow(e,4)-(288000-115200*eta)*e*e
			+57600-23040*eta)+2*(-(3240-8520*eta-3960*eta*eta)*pow(ecosu,5)+(12960-34080*eta-15840*eta*eta)*pow(ecosu,4)+((6580-17040*eta-7920*eta*eta)*e*e-12582+24070*eta
			+43770*eta*eta)*pow(ecosu,3)+(-(35820-101340*eta-21060*eta*eta)*e*e+6846+4530*eta-78930*eta*eta)*ecosu*ecosu+((28455-74415*eta-33825*eta*eta)*e*e+12159-64695*eta
			+85575*eta*eta)*ecosu+(20655-68535*eta+21735*eta*eta)*pow(e,4)-(40425-127185*eta+22785*eta*eta)*e*e+4512-6880*eta-16800*eta*eta)*C*C+15*(1-5*eta
			+5*eta*eta)*(120*pow(ecosu,5)+480*pow(ecosu,4)+(240*e*e-806)*pow(ecosu,3)-(900*e*e-918)*ecosu*ecosu+(955*e*e-613)*ecosu-45*pow(e,4)-205*e*e+96)*pow(C,4));
	return term;
}
static double  X_2_C_4_C_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S*S/24*esinu*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-1./2.)*(((120-344*eta-72*eta*eta)*e*e-72+152*eta+216*eta*eta)*pow(ecosu,3)+(-(480-1400*eta-192*eta*eta)*e*e
			+324-776*eta-660*eta*eta)*ecosu*ecosu+(-(240-688*eta-144*eta*eta)*pow(e,4)+(518-1126*eta-1442*eta*eta)*e*e-110-234*eta+1802*eta*eta)*ecosu
			+(831-2425*eta-147*eta*eta)*pow(e,4)-(1340-3544*eta-1328*eta*eta)*e*e+449-879*eta-1361*eta*eta)+Bpow(1-e*e,1./2.)*(1-5*eta+5*eta*eta)*(-48*pow(ecosu,3)
			+162*ecosu*ecosu+(96*e*e-250)*ecosu-201*e*e+241)*C*C);
	return term;
}
static double  X_2_C_4_S_4(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S*S/720*Bpow(1-ecosu,-5)*(((900-2580*eta-540*eta*eta)*pow(ecosu,5)-(3600-10320*eta-2160*eta*eta)*pow(ecosu,4)+(-(7200-20640*eta-4320*eta*eta)*e*e
			+4830-7240*eta-24120*eta*eta)*pow(ecosu,3)+((27990-82350*eta-7290*eta*eta)*e*e-12120+21390*eta+47610*eta*eta)*ecosu*ecosu+((7200-20640*eta
			-4320*eta*eta)*pow(e,4)-(8430-12105*eta-42675*eta*eta)*e*e-13500+61875*eta-65535*eta*eta)*ecosu-(24930-74325*eta+3465*eta*eta)*pow(e,4)
			+(23100-57765*eta-24135*eta*eta)*e*e+5760-30080*eta+32640*eta*eta)+3*(1-5*eta+5*eta*eta)*(120*pow(ecosu,5)-480*pow(ecosu,4)+(1526-960*e*e)*pow(ecosu,3)
			+(2700*e*e-2718)*ecosu*ecosu+(960*pow(e,4)-3235*e*e+1933)*ecosu-2115*pow(e,4)+3805*e*e-1536)*C*C);
	return term;
}
static double  X_2_C_6_C_6(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*pow(S,4)/96*(1-5*eta+5*eta*eta)*Bpow(1-e*e,1./2.)*esinu*Bpow(1-ecosu,-6)*(-72*pow(ecosu,4)+312*pow(ecosu,3)+(384*e*e-844)*ecosu*ecosu
			-(1053*e*e-1325)*ecosu-384*pow(e,4)+1437*e*e-1105);
	return term;
}
static double  X_2_C_6_S_6(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*pow(S,4)/960*(1-5*eta+5*eta*eta)*Bpow(1-ecosu,-6)*(120*pow(ecosu,6)-600*pow(ecosu,5)-(2160*e*e-3206)*pow(ecosu,4)+(7860*e*e-8444)*pow(ecosu,3)
			+(5760*pow(e,4)-19135*e*e+13051)*ecosu*ecosu-(11475*pow(e,4)-23240*e*e+11269)*ecosu-3840*pow(e,6)+17235*pow(e,4)-21325*e*e+7776);
	return term;
}
static double  X_2(double e, double u, double iota, double eta){
	double ecosu = e*cos(u);
	double esinu = e*sin(u);
	double C = cos(iota);
	double S = sin(iota);
	double term;

	term = C*S*S/24*esinu*Bpow(1-ecosu,-5)*(Bpow(1-e*e,-1./2.)*((-(36-116*eta+12*eta*eta)*e*e+24-68*eta-24*eta*eta)*ecosu*ecosu+((174-538*eta+10*eta*eta)*e*e
			-150+442*eta+62*eta*eta)*ecosu-(153-459*eta-21*eta*eta)*pow(e,4)+(168-496*eta-40*eta*eta)*e*e-27+85*eta-17*eta*eta)+Bpow(1-e*e,1./2.)*(1
			+5*eta*eta-5*eta)*(5*ecosu*ecosu-22*ecosu+15*e*e+1)*C*C);
	return term;
}
//below are the different PN order H plus and cross (note the PN expansion parameter is not included. These follow eqns 2.24 and 2.25 of gopu and Iyer (2002))

static double PN2_plus_gauge_corr(double e, double u, double w, \
double lam, double iota, double eta)
{
return (e*(1 + 17*eta)*Bpow(-1 + \
e*cos(u),-3)*(Bpow(sin(iota),2)*cos(u)*(-1 + e*cos(u)) - (1 + \
Bpow(cos(iota),2))*cos(2*lam)*(-4*e + e*Bpow(cos(u),2) + \
3*cos(u))*cos(2*w) + (1 + Bpow(cos(iota),2))*(-4*e + e*Bpow(cos(u),2) \
+ 3*cos(u))*sin(2*lam)*sin(2*w) - 2*Bpow(1 - Bpow(e,2),-0.5)*(1 + \
Bpow(cos(iota),2))*(1 - 2*Bpow(e,2) + \
e*cos(u))*sin(u)*(cos(2*w)*sin(2*lam) + cos(2*lam)*sin(2*w))))/4.;
}

static double PN2_cross_gauge_corr(double e, double u, double w, \
double lam, double iota, double eta)
{
return -(e*(1 + 17*eta)*Bpow(1 - Bpow(e,2),-0.5)*Bpow(-1 + \
e*cos(u),-3)*cos(iota)*(cos(2*w)*(-4*e*Bpow(1 - \
Bpow(e,2),0.5)*sin(2*lam) + e*Bpow(1 - \
Bpow(e,2),0.5)*Bpow(cos(u),2)*sin(2*lam) + 3*Bpow(1 - \
Bpow(e,2),0.5)*cos(u)*sin(2*lam) - 2*cos(2*lam)*sin(u) + \
4*Bpow(e,2)*cos(2*lam)*sin(u) - 2*e*cos(2*lam)*cos(u)*sin(u)) + \
Bpow(1 - Bpow(e,2),0.5)*cos(2*lam)*(-4*e + e*Bpow(cos(u),2) + \
3*cos(u))*sin(2*w) + 2*(1 - 2*Bpow(e,2) + \
e*cos(u))*sin(2*lam)*sin(u)*sin(2*w)))/2.;
}

static double  H_P_0(double e, double u, double w, double lam, double iota){
	double P0C2C2 =  P_0_C_2_C_2(e, u, iota);
	double P0C2S2 =  P_0_C_2_S_2(e, u, iota);
	double P0 = P_0(e, u, iota);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double term;

	term = (P0C2C2*c2w+P0C2S2*s2w)*c2lam+(P0C2S2*c2w-P0C2C2*s2w)*s2lam+P0;
	return term;
}
static double  H_P_05(double e, double u, double w, double lam, double iota, double m1, double m2){
	double P05C1C1 = P_05_C_1_C_1(e, u, iota, m1, m2);
	double P05C1S1 = P_05_C_1_S_1(e, u, iota, m1, m2);
	double P05C3C3 = P_05_C_3_C_3(e, u, iota, m1, m2);
	double P05C3S3 = P_05_C_3_S_3(e, u, iota, m1, m2);
	double cw = cos(w);
	double sw = sin(w);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c3w = cos(3*w);
	double s3w = sin(3*w);
	double clam = cos(lam);
	double slam = sin(lam);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c3lam = cos(3*lam);
	double s3lam = sin(3*lam);
	double term;

	term = (P05C1C1*cw+P05C1S1*sw)*clam+(P05C1S1*cw-P05C1C1*sw)*slam+(P05C3C3*c3w+P05C3S3*s3w)*c3lam+(P05C3S3*c3w-P05C3C3*s3w)*s3lam;
	return term;
}
static double  H_P_1(double e, double u, double w, double lam, double iota, double eta){
	double P1C2C2 = P_1_C_2_C_2(e, u, iota, eta);
	double P1C2S2 = P_1_C_2_S_2(e, u, iota, eta);
	double P1C4C4 = P_1_C_4_C_4(e, u, iota, eta);
	double P1C4S4 = P_1_C_4_S_4(e, u, iota, eta);
	double P1 = P_1(e, u, iota, eta);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c4w = cos(4*w);
	double s4w = sin(4*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c4lam = cos(4*lam);
	double s4lam = sin(4*lam);
	double term;

	term = (P1C2C2*c2w+P1C2S2*s2w)*c2lam+(P1C2S2*c2w-P1C2C2*s2w)*s2lam+(P1C4C4*c4w+P1C4S4*s4w)*c4lam+(P1C4S4*c4w-P1C4C4*s4w)*s4lam+P1;
	return term;
}
static double  H_P_15(double e, double u, double w, double lam, double iota, double m1, double m2, double eta){
	double P15C1C1 = P_15_C_1_C_1(e, u, iota, m1, m2, eta);
	double P15C1S1 = P_15_C_1_S_1(e, u, iota, m1, m2, eta);
	double P15C3C3 = P_15_C_3_C_3(e, u, iota, m1, m2, eta);
	double P15C3S3 = P_15_C_3_S_3(e, u, iota, m1, m2, eta);
	double P15C5C5 = P_15_C_5_C_5(e, u, iota, m1, m2, eta);
	double P15C5S5 = P_15_C_5_S_5(e, u, iota, m1, m2, eta);
	double cw = cos(w);
	double sw = sin(w);
	double c3w = cos(3*w);
	double s3w = sin(3*w);
	double c5w = cos(5*w);
	double s5w = sin(5*w);
	double clam = cos(lam);
	double slam = sin(lam);
	double c3lam = cos(3*lam);
	double s3lam = sin(3*lam);
	double c5lam = cos(5*lam);
	double s5lam = sin(5*lam);
	double term;

	term = (P15C1C1*cw+P15C1S1*sw)*clam+(P15C1S1*cw-P15C1C1*sw)*slam+(P15C3C3*c3w+P15C3S3*s3w)*c3lam+(P15C3S3*c3w-P15C3C3*s3w)*s3lam
			+(P15C5C5*c5w-P15C3S3*s5w)*c5lam+(P15C5S5*c5w-P15C5C5*s5w)*s5lam;
	return term;
}
static double  H_P_2(double e, double u, double w, double lam, double iota, double eta){
	double P2C2C2 = P_2_C_2_C_2(e, u, iota, eta);
	double P2C2S2 = P_2_C_2_S_2(e, u, iota, eta);
	double P2C4C4 = P_2_C_4_C_4(e, u, iota, eta);
	double P2C4S4 = P_2_C_4_S_4(e, u, iota, eta);
	double P2C6C6 = P_2_C_6_C_6(e, u, iota, eta);
	double P2C6S6 = P_2_C_6_S_6(e, u, iota, eta);
	double P2 = P_2(e, u, iota, eta);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c4w = cos(4*w);
	double s4w = sin(4*w);
	double c6w = cos(6*w);
	double s6w = sin(6*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c4lam = cos(4*lam);
	double s4lam = sin(4*lam);
	double c6lam = cos(6*lam);
	double s6lam = sin(6*lam);
	double term;
	double gaugecorr = PN2_plus_gauge_corr(e, u, w, lam, iota, eta);

	term = (P2C2C2*c2w+P2C2S2*s2w)*c2lam+(P2C2S2*c2w-P2C2C2*s2w)*s2lam+(P2C4C4*c4w+P2C4S4*s4w)*c4lam+(P2C4S4*c4w-P2C4C4*s4w)*s4lam+(P2C6C6*c6w
			+P2C6S6*s6w)*c6lam+(P2C6S6*c6w-P2C6C6*s6w)*s6lam+P2;
	return term + gaugecorr;
}
static double  H_X_0(double e, double u, double w, double lam, double iota){
	double X0C2C2 = X_0_C_2_C_2(e, u, iota);
	double X0C2S2 = X_0_C_2_S_2(e, u, iota);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double term;

	term = (X0C2C2*c2w+X0C2S2*s2w)*c2lam+(X0C2S2*c2w-X0C2C2*s2w)*s2lam;
	return term;
}
static double  H_X_05(double e, double u, double w, double lam, double iota, double m1, double m2){
	double X05C1C1 = X_05_C_1_C_1(e, u, iota, m1, m2);
	double X05C1S1 = X_05_C_1_S_1(e, u, iota, m1, m2);
	double X05C3C3 = X_05_C_3_C_3(e, u, iota, m1, m2);
	double X05C3S3 = X_05_C_3_S_3(e, u, iota, m1, m2);
	double cw = cos(w);
	double sw = sin(w);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c3w = cos(3*w);
	double s3w = sin(3*w);
	double clam = cos(lam);
	double slam = sin(lam);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c3lam = cos(3*lam);
	double s3lam = sin(3*lam);
	double term;

	term = (X05C1C1*cw+X05C1S1*sw)*clam+(X05C1S1*cw-X05C1C1*sw)*slam+(X05C3C3*c3w+X05C3S3*s3w)*c3lam+(X05C3S3*c3w-X05C3C3*s3w)*s3lam;
	return term;
}
static double  H_X_1(double e, double u, double w, double lam, double iota, double eta){
	double X1C2C2 = X_1_C_2_C_2(e ,u, iota, eta);
	double X1C2S2 = X_1_C_2_S_2(e, u, iota, eta);
	double X1C4C4 = X_1_C_4_C_4(e, u, iota, eta);
	double X1C4S4 = X_1_C_4_S_4(e, u, iota, eta);
	double X1 = X_1(e, u, iota, eta);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c4w = cos(4*w);
	double s4w = sin(4*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c4lam = cos(4*lam);
	double s4lam = sin(4*lam);
	double term;

	term = (X1C2C2*c2w+X1C2S2*s2w)*c2lam+(X1C2S2*c2w-X1C2C2*s2w)*s2lam+(X1C4C4*c4w+X1C4S4*s4w)*c4lam+(X1C4S4*c4w-X1C4C4*s4w)*s4lam+X1;
	return term;
}
static double  H_X_15(double e, double u, double w, double lam, double iota, double m1, double m2, double eta){
	double X15C1C1 = X_15_C_1_C_1(e, u, iota, m1, m2, eta);
	double X15C1S1 = X_15_C_1_S_1(e, u, iota, m1, m2, eta);
	double X15C3C3 = X_15_C_3_C_3(e, u, iota, m1, m2, eta);
	double X15C3S3 = X_15_C_3_S_3(e, u, iota, m1, m2, eta);
	double X15C5C5 = X_15_C_5_C_5(e, u, iota, m1, m2, eta);
	double X15C5S5 = X_15_C_5_S_5(e, u, iota, m1, m2, eta);
	double cw = cos(w);
	double sw = sin(w);
	double c3w = cos(3*w);
	double s3w = sin(3*w);
	double c5w = cos(5*w);
	double s5w = sin(5*w);
	double clam = cos(lam);
	double slam = sin(lam);
	double c3lam = cos(3*lam);
	double s3lam = sin(3*lam);
	double c5lam = cos(5*lam);
	double s5lam = sin(5*lam);
	double term;

	term = (X15C1C1*cw+X15C1S1*sw)*clam+(X15C1S1*cw-X15C1C1*sw)*slam+(X15C3C3*c3w+X15C3S3*s3w)*c3lam+(X15C3S3*c3w-X15C3C3*s3w)*s3lam
			+(X15C5C5*c5w+X15C5S5*s5w)*c5lam+(X15C5S5*c5w-X15C5C5*s5w)*s5lam;
	return term;
}
static double  H_X_2(double e, double u, double w, double lam, double iota, double eta){
	double X2C2C2 = X_2_C_2_C_2(e, u, iota, eta);
	double X2C2S2 = X_2_C_2_S_2(e, u, iota, eta);
	double X2C4C4 = X_2_C_4_C_4(e, u, iota, eta);
	double X2C4S4 = X_2_C_4_S_4(e, u, iota, eta);
	double X2C6C6 = X_2_C_6_C_6(e, u, iota, eta);
	double X2C6S6 = X_2_C_6_S_6(e, u, iota, eta);
	double X2 = X_2(e, u, iota, eta);
	double c2w = cos(2*w);
	double s2w = sin(2*w);
	double c4w = cos(4*w);
	double s4w = sin(4*w);
	double c6w = cos(6*w);
	double s6w = sin(6*w);
	double c2lam = cos(2*lam);
	double s2lam = sin(2*lam);
	double c4lam = cos(4*lam);
	double s4lam = sin(4*lam);
	double c6lam = cos(6*lam);
	double s6lam = sin(6*lam);
	double term;
	double gaugecorr = PN2_cross_gauge_corr(e, u, w, lam, iota, eta);

	term = (X2C2C2*c2w+X2C2S2*s2w)*c2lam+(X2C2S2*c2w-X2C2C2*s2w)*s2lam+(X2C4C4*c4w+X2C4S4*s4w)*c4lam+(X2C4S4*c4w-X2C4C4*s4w)*s4lam
			+(X2C6C6*c6w+X2C6S6*s6w)*c6lam+(X2C6S6*c6w-X2C6C6*s6w)*s6lam+X2;
	return term + gaugecorr;
}


#endif /* TD_PN_AMPS_HPP_ */
