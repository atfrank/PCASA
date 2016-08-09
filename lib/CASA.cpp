//Aaron T. Frank

/*
This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

//Code generated using: awk '{print "{\""$1":"$2":"$3"\","$4"},"}' larmorD_both.dat | tr '\n' ' '

#include "CASA.hpp"
#include "Molecule.hpp"
#include "Misc.hpp"

#include <fstream>
#include <stdlib.h>
#include <algorithm>


CASA::CASA (Molecule *mol, const std::string fparmfile, const std::string freffile, const std::string fsasa){
    this->initializeAminoAcids();        
    this->initializeHist();

    /* load SASA info from file */
    if (fsasa.length() > 0){
      this->loadSASAFile(fsasa);
    }
    /* load parameters */
    if (fparmfile.length() > 0){
      this->loadParmFile(fparmfile);
    } 
    else {
    	this->initializeSASACoef();
    }
    /* load reference SASA */
    if (freffile.length() > 0){
      this->loadRefFile(freffile);
    }
    else {
    	this->initializeRefSASA();
    }
}

void CASA::initializeAminoAcids(){
  this->aminoAcids.push_back("ALA");
  this->aminoAcids.push_back("ARG");
  this->aminoAcids.push_back("ASN");
  this->aminoAcids.push_back("ASP");
  this->aminoAcids.push_back("CYS");
  this->aminoAcids.push_back("GLN");
  this->aminoAcids.push_back("GLU");
  this->aminoAcids.push_back("GLY");
  this->aminoAcids.push_back("HIS");
  this->aminoAcids.push_back("ILE");
  this->aminoAcids.push_back("LEU");
  this->aminoAcids.push_back("LYS");
  this->aminoAcids.push_back("MET");
  this->aminoAcids.push_back("PHE");
  this->aminoAcids.push_back("PRO");
  this->aminoAcids.push_back("SER");
  this->aminoAcids.push_back("THR");
  this->aminoAcids.push_back("TRP");
  this->aminoAcids.push_back("TYR");
  this->aminoAcids.push_back("VAL");  
  std::vector<std::string>::iterator it=std::unique(this->aminoAcids.begin(), this->aminoAcids.end());
  this->aminoAcids.resize(std::distance(this->aminoAcids.begin(),it));
}

void CASA::initializeRefSASA(){
  this->refSASA.insert(std::pair<std::string,double>("ALA",137.50496643251));
  this->refSASA.insert(std::pair<std::string,double>("ARG",236.902235217962));
  this->refSASA.insert(std::pair<std::string,double>("ASN",184.134864518706));
  this->refSASA.insert(std::pair<std::string,double>("ASP",186.734811264804));
  this->refSASA.insert(std::pair<std::string,double>("CYS",153.447234418197));
  this->refSASA.insert(std::pair<std::string,double>("GLN",220.199142372231));
  this->refSASA.insert(std::pair<std::string,double>("GLU",216.957744464548));
  this->refSASA.insert(std::pair<std::string,double>("GLY",101.908174205877));
  this->refSASA.insert(std::pair<std::string,double>("HIS",183.774948105238));
  this->refSASA.insert(std::pair<std::string,double>("ILE",211.179667274773));
  this->refSASA.insert(std::pair<std::string,double>("LEU",189.262346285463));
  this->refSASA.insert(std::pair<std::string,double>("LYS",206.1419097756));
  this->refSASA.insert(std::pair<std::string,double>("MET",208.695656914473));
  this->refSASA.insert(std::pair<std::string,double>("PHE",261.102510143459));
  this->refSASA.insert(std::pair<std::string,double>("PRO",178.325572209738));
  this->refSASA.insert(std::pair<std::string,double>("SER",150.044108409222));
  this->refSASA.insert(std::pair<std::string,double>("THR",173.029721961006));
  this->refSASA.insert(std::pair<std::string,double>("TRP",270.4835515474));
  this->refSASA.insert(std::pair<std::string,double>("TYR",265.456513268642));
  this->refSASA.insert(std::pair<std::string,double>("VAL",185.401997195671));
}

void CASA::initializeSASACoef(){
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:ARG", -83.8207772727775));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:GLY", -42.3643376774708));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:ALA", -57.466627380429));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:SER", -62.5417094241778));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:CYS", -60.6472376416871));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:PRO", -75.1440566921417));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:THR", -72.0795478271455));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:ASP", -73.8841310987875));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:VAL", -73.0255824038019));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:ASN", -68.4338959150843));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:GLU", -78.416772320484));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:ILE", -81.0121707552322));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:LEU", -80.5315806988284));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:GLN", -77.624563954462));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:HIS", -71.2619264411202));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:MET", -71.8203519265687));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:LYS", -79.5833667736523));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:PHE", -84.959899934399));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:TYR", -87.2065731803885));
	this->sasaCoef.insert(std::pair<std::string,double>("ARG:TRP", -87.0012943979153));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:ARG", -42.0261201004127));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:GLY", -17.6998949959275));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:ALA", -22.2473870887928));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:SER", -25.1841957696804));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:CYS", -26.0965565937699));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:PRO", -35.4777967551648));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:THR", -29.4293365655117));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:ASP", -30.7121495549564));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:VAL", -29.7391636952366));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:ASN", -29.8018245104223));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:GLU", -35.169183715146));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:ILE", -33.3849879909354));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:LEU", -33.6463492374973));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:GLN", -34.1726619079028));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:HIS", -36.0748905171377));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:MET", -36.0968694781901));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:LYS", -38.9443518052504));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:PHE", -35.9320925119826));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:TYR", -37.8841728508113));
	this->sasaCoef.insert(std::pair<std::string,double>("GLY:TRP", -40.3410346578774));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:ARG", -53.6960627222511));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:GLY", -21.8545974175632));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:ALA", -33.0825688139907));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:SER", -36.5059419317034));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:CYS", -34.7500843688118));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:PRO", -42.8468928278172));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:THR", -39.5482207583115));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:ASP", -41.8825246622213));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:VAL", -41.8161216613512));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:ASN", -40.2881532517147));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:GLU", -47.6184489421926));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:ILE", -47.307691277721));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:LEU", -46.4196639357636));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:GLN", -49.7794678346143));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:HIS", -45.0884537041172));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:MET", -46.0730215544936));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:LYS", -51.2833020945776));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:PHE", -49.3057428745254));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:TYR", -48.5905784222719));
	this->sasaCoef.insert(std::pair<std::string,double>("ALA:TRP", -52.5743947010102));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:ARG", -59.8017232127831));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:GLY", -25.5954690237441));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:ALA", -35.0447567701382));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:SER", -37.6327766140352));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:CYS", -34.3426091524195));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:PRO", -49.1789425577898));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:THR", -46.9215102685333));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:ASP", -47.4563887555813));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:VAL", -44.6694747035461));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:ASN", -40.9890409901454));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:GLU", -52.9029876268529));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:ILE", -50.672980534798));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:LEU", -51.8651344702657));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:GLN", -51.575721945785));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:HIS", -46.1363753214778));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:MET", -51.3937786867948));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:LYS", -54.5063367407581));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:PHE", -54.0362412451941));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:TYR", -54.7714189944317));
	this->sasaCoef.insert(std::pair<std::string,double>("SER:TRP", -54.9732325941883));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:ARG", -58.4382142698828));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:GLY", -27.2235575336866));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:ALA", -38.4278497819029));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:SER", -40.0065906124384));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:CYS", -38.3657468024512));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:PRO", -51.7663618493384));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:THR", -42.9110152107167));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:ASP", -44.6807918423797));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:VAL", -46.3715274123085));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:ASN", -45.7571769517819));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:GLU", -51.2782075338388));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:ILE", -52.5123941310801));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:LEU", -52.3623985602548));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:GLN", -47.4622035140449));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:HIS", -49.8556239673118));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:MET", -42.6149212555577));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:LYS", -55.6186823226335));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:PHE", -52.4072056922846));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:TYR", -47.5464148312276));
	this->sasaCoef.insert(std::pair<std::string,double>("CYS:TRP", -48.1968248654214));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:ARG", -72.0462977658937));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:GLY", -26.5887601249948));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:ALA", -41.8794531465695));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:SER", -43.8269886894937));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:CYS", -37.9919945688546));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:PRO", -55.2497435344688));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:THR", -53.6892451715948));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:ASP", -53.0105705617222));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:VAL", -56.6073062121394));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:ASN", -48.6916388085815));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:GLU", -57.6771696554508));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:ILE", -66.27176649793));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:LEU", -61.0925839966704));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:GLN", -58.6349084610492));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:HIS", -55.246281191655));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:MET", -50.9684251638297));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:LYS", -61.8609972821172));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:PHE", -64.0021063989349));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:TYR", -57.5831150556712));
	this->sasaCoef.insert(std::pair<std::string,double>("PRO:TRP", -65.8400863952117));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:ARG", -66.6757378621554));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:GLY", -31.7127561075626));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:ALA", -41.4947771057919));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:SER", -41.3904818334398));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:CYS", -37.9860252376184));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:PRO", -56.4877760170235));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:THR", -46.9312035125275));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:ASP", -56.1154985941805));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:VAL", -50.4371564733685));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:ASN", -51.3919755637858));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:GLU", -58.5057446802779));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:ILE", -56.4628589772851));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:LEU", -58.6239529730061));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:GLN", -61.4470403583784));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:HIS", -57.8921070547608));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:MET", -60.6646945137906));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:LYS", -60.73580311225));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:PHE", -62.7537983897095));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:TYR", -61.6960309578404));
	this->sasaCoef.insert(std::pair<std::string,double>("THR:TRP", -64.9509343888304));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:ARG", -70.6773782369662));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:GLY", -34.5263041508645));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:ALA", -47.9386812841558));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:SER", -47.2665505242678));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:CYS", -44.1302164223373));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:PRO", -61.535620829317));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:THR", -57.0457027164369));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:ASP", -59.2360196668507));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:VAL", -55.7287876965179));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:ASN", -59.9238956997147));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:GLU", -62.0133799248102));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:ILE", -62.2709567580556));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:LEU", -62.4428791120909));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:GLN", -60.2737947222094));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:HIS", -58.6397120050227));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:MET", -56.3225518997233));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:LYS", -61.5153190062843));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:PHE", -66.0429936110929));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:TYR", -66.8159420310155));
	this->sasaCoef.insert(std::pair<std::string,double>("ASP:TRP", -70.0396266776663));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:ARG", -71.8557697407599));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:GLY", -32.3534386335658));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:ALA", -44.0716175762182));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:SER", -49.5542063020519));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:CYS", -47.2584079841836));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:PRO", -56.8876261462732));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:THR", -56.6367644588696));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:ASP", -59.810293044309));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:VAL", -57.4458341972382));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:ASN", -57.3710166788525));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:GLU", -64.5961377968235));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:ILE", -60.7261243150928));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:LEU", -63.6525500012758));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:GLN", -66.4431942041292));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:HIS", -59.2530371265182));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:MET", -61.1484189063272));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:LYS", -64.8590870960905));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:PHE", -65.2782975386906));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:TYR", -66.6066534179694));
	this->sasaCoef.insert(std::pair<std::string,double>("VAL:TRP", -69.0745701540431));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:ARG", -69.4531847780047));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:GLY", -37.2714326736));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:ALA", -48.2863793966305));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:SER", -48.3933687555913));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:CYS", -38.8940155700208));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:PRO", -57.2842337127655));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:THR", -52.5694439255301));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:ASP", -60.2489794760674));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:VAL", -55.9722716600862));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:ASN", -53.1768065488211));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:GLU", -64.6958985500376));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:ILE", -59.2838699191939));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:LEU", -63.1062376485673));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:GLN", -62.7868221359032));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:HIS", -59.3931105540753));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:MET", -54.618233144183));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:LYS", -61.9997684239135));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:PHE", -64.6886500206045));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:TYR", -65.331946855965));
	this->sasaCoef.insert(std::pair<std::string,double>("ASN:TRP", -68.8678177371111));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:ARG", -78.9692957776615));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:GLY", -44.6461752548837));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:ALA", -52.8351227319165));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:SER", -58.3161777549291));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:CYS", -53.442997397828));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:PRO", -69.6214511025424));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:THR", -67.2869730032662));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:ASP", -64.2373671447865));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:VAL", -68.1168638350652));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:ASN", -61.797303402986));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:GLU", -67.2248068495065));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:ILE", -73.799591825698));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:LEU", -71.5388415724329));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:GLN", -69.1582388809183));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:HIS", -67.6682996289317));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:MET", -69.0620888370368));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:LYS", -69.4192966635131));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:PHE", -80.0088441181985));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:TYR", -76.7572220177097));
	this->sasaCoef.insert(std::pair<std::string,double>("GLU:TRP", -79.2849104101016));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:ARG", -81.4777140568317));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:GLY", -37.4649412612081));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:ALA", -50.872852394874));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:SER", -58.030651948003));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:CYS", -45.9031477033497));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:PRO", -67.3507411717135));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:THR", -64.7840174861555));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:ASP", -64.892581595724));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:VAL", -66.2129513937875));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:ASN", -63.38871510859));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:GLU", -69.9852660773031));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:ILE", -73.0525149142362));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:LEU", -73.5157655343603));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:GLN", -73.0787527646084));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:HIS", -67.1925522258214));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:MET", -64.6836867085342));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:LYS", -74.3930682872983));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:PHE", -77.4691238757378));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:TYR", -75.3939393559114));
	this->sasaCoef.insert(std::pair<std::string,double>("ILE:TRP", -76.9483121540776));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:ARG", -72.5844768320253));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:GLY", -33.9550266534096));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:ALA", -46.9035198960781));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:SER", -49.340068744134));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:CYS", -45.4286771523816));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:PRO", -60.824809461367));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:THR", -53.7687962489564));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:ASP", -58.4932045252552));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:VAL", -58.2009398376905));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:ASN", -57.6024693521829));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:GLU", -69.3409168363889));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:ILE", -62.3785369370504));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:LEU", -64.4832711266644));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:GLN", -67.013440497649));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:HIS", -60.4887712824995));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:MET", -61.0919213546262));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:LYS", -68.5119344388417));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:PHE", -71.6466274456572));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:TYR", -67.8944756710546));
	this->sasaCoef.insert(std::pair<std::string,double>("LEU:TRP", -69.5597454048553));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:ARG", -77.5156163258467));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:GLY", -46.1524609383607));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:ALA", -52.5636932607279));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:SER", -58.8843597471767));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:CYS", -55.0318398496677));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:PRO", -71.0413283471267));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:THR", -66.2944607622895));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:ASP", -66.2196048749543));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:VAL", -69.2642563391464));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:ASN", -66.4844378300654));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:GLU", -68.4135847365013));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:ILE", -72.5817870475189));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:LEU", -71.8291838731274));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:GLN", -67.596698992711));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:HIS", -66.8901091026169));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:MET", -66.5257156543246));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:LYS", -72.9529787059551));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:PHE", -75.2203431330683));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:TYR", -72.8170668327023));
	this->sasaCoef.insert(std::pair<std::string,double>("GLN:TRP", -79.9826539534809));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:ARG", -69.1223814025902));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:GLY", -29.6496296906568));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:ALA", -46.9449211862017));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:SER", -47.294421804825));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:CYS", -40.4181606995592));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:PRO", -56.3120944791406));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:THR", -51.8766436259439));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:ASP", -53.4026211423083));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:VAL", -55.1402070018892));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:ASN", -51.0799703158856));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:GLU", -60.2351362437989));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:ILE", -59.7033852495159));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:LEU", -61.7495435629866));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:GLN", -66.3915662560925));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:HIS", -50.9396829245052));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:MET", -56.9920566978722));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:LYS", -64.9556407312554));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:PHE", -68.974944283962));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:TYR", -60.9883879251745));
	this->sasaCoef.insert(std::pair<std::string,double>("HIS:TRP", -62.419537237404));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:ARG", -80.6208465499368));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:GLY", -42.5024870804051));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:ALA", -52.3791691920921));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:SER", -54.7205482714701));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:CYS", -47.0397313722078));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:PRO", -62.4561651815219));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:THR", -64.9033142805557));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:ASP", -66.3180158592913));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:VAL", -63.4396135186527));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:ASN", -60.1476397898014));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:GLU", -71.8396466450588));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:ILE", -68.2305829127021));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:LEU", -71.0672072692088));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:GLN", -65.7409926709568));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:HIS", -58.8441603386642));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:MET", -61.0215405428503));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:LYS", -67.0451588284521));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:PHE", -75.4320233394029));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:TYR", -73.6717869032781));
	this->sasaCoef.insert(std::pair<std::string,double>("MET:TRP", -66.7117466322911));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:ARG", -77.6253778833452));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:GLY", -38.0427543999273));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:ALA", -48.3387210020821));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:SER", -49.9576488329855));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:CYS", -47.548853547874));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:PRO", -63.3080137659962));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:THR", -61.8981454119626));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:ASP", -63.2595828067634));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:VAL", -66.5998137768474));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:ASN", -64.7473051402651));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:GLU", -66.9391609217043));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:ILE", -73.1637036718286));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:LEU", -73.1442870109027));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:GLN", -65.6236194317787));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:HIS", -63.2638773701353));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:MET", -64.4397564802315));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:LYS", -67.0906045397324));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:PHE", -73.5217752073495));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:TYR", -69.7445544286352));
	this->sasaCoef.insert(std::pair<std::string,double>("LYS:TRP", -80.6728587376575));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:ARG", -92.6730083465509));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:GLY", -47.4476465004194));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:ALA", -65.8493162066439));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:SER", -71.7408364861269));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:CYS", -60.1548057936917));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:PRO", -81.7956337452181));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:THR", -77.8532378560295));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:ASP", -78.7067926639052));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:VAL", -82.3555876646174));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:ASN", -75.8125690258605));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:GLU", -86.7326694845289));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:ILE", -91.6177025095725));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:LEU", -92.1786439746122));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:GLN", -81.2600900505726));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:HIS", -75.7745255227411));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:MET", -79.2748759423615));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:LYS", -88.8973067129848));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:PHE", -95.6536511713865));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:TYR", -92.5523163160029));
	this->sasaCoef.insert(std::pair<std::string,double>("PHE:TRP", -91.0037368821828));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:ARG", -90.088666913567));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:GLY", -50.6025499936813));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:ALA", -64.4912565750068));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:SER", -67.8524743463385));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:CYS", -55.7217179293744));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:PRO", -77.7694890307798));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:THR", -77.8216597026619));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:ASP", -77.7946030020674));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:VAL", -84.2438798502286));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:ASN", -77.8010677051441));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:GLU", -84.7308569102312));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:ILE", -90.9937974845723));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:LEU", -90.1553833445557));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:GLN", -76.6296935764277));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:HIS", -75.2097108881859));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:MET", -72.5395732961266));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:LYS", -81.8124625790141));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:PHE", -95.8969624552547));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:TYR", -90.454588608311));
	this->sasaCoef.insert(std::pair<std::string,double>("TYR:TRP", -79.960404694134));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:ARG", -80.4143132629891));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:GLY", -55.4199519048124));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:ALA", -69.6189837929677));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:SER", -68.1595795971853));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:CYS", -47.3644870927432));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:PRO", -64.4145759079018));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:THR", -76.423676165284));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:ASP", -71.5171312087642));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:VAL", -85.1380429521876));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:ASN", -67.690939431147));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:GLU", -74.7671003133502));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:ILE", -90.0974572199688));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:LEU", -91.2931611893283));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:GLN", -64.7678475961554));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:HIS", -57.5241722265667));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:MET", -58.4961804503981));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:LYS", -76.311359460891));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:PHE", -96.0514761698289));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:TYR", -86.2232150707108));
	this->sasaCoef.insert(std::pair<std::string,double>("TRP:TRP", -66.5390468648662));
}

void CASA::initializeHist(){
  for (unsigned int i=0; i< this->aminoAcids.size(); i++){
    this->hist.insert(std::pair<std::string,double>(this->aminoAcids.at(i),0.0));
  }
}

void CASA::setHist(const std::string &key, const double &val){
    if (this->hist.find (key) != this->hist.end()){
      this->hist.at(key) = val;
    } else {
      this->hist.insert(std::pair<std::string, double>(key, val));
    }
}

void CASA::accumHist(const std::string &key, const double &val){
    if (this->hist.find (key) != this->hist.end()){
      this->hist.at(key) += val;
    } else {
      this->hist.insert(std::pair<std::string, double>(key, val));
    }
    
}

double CASA::predictHist(const std::string &resname){
  double predSASA = 0.0;
  std::string res = "";
  std::string key = "";
  for (unsigned int i=0; i< this->aminoAcids.size(); i++){
    res = this->aminoAcids.at(i);
    key = resname+":"+res;
    predSASA += this->sasaCoef.at(key)*this->hist.at(res);
  }
  return(predSASA);
}

void CASA::printHist(){
  for (unsigned int i=0; i< this->aminoAcids.size(); i++){
    std::cout << this->hist.at(this->aminoAcids.at(i)) << " ";
  }
}

void CASA::printHistHeader(){
  for (unsigned int i=0; i< this->aminoAcids.size(); i++){
    std::cout << this->aminoAcids.at(i) << " ";
  }
}

void CASA::clearHist(){
    this->hist.clear();
}

double CASA::getRefSASA(const std::string &key){
    if (this->refSASA.find (key) == this->refSASA.end()){
        return 0.0;
    } else {
        return (this->refSASA.at(key));
    }
}

std::string CASA::printRefSASAVec(const std::string &key){    
    std::stringstream s;
    s.str("");
    for (unsigned int i=0; i< this->aminoAcids.size(); i++){
    	if (this->aminoAcids.at(i) == key){
    		 s << this->refSASA.at(key) << " ";
    	}
    	else {
    		s << "0.0 ";
    	}
    }
    return(s.str());
}

std::string CASA::printRefSASAVecHeader(){    
	std::stringstream s;
	s.str("");
	for (unsigned int i=0; i< this->aminoAcids.size(); i++){
		s << "RefSASA_" << aminoAcids.at(i) << " ";
	}
	return(s.str());
}

double CASA::getSASACoef(const std::string &key){
    if (this->sasaCoef.find (key) == this->sasaCoef.end()){
        return 0.0;
    } else {
        return (this->sasaCoef.at(key));
    }
}

double CASA::getHist(const std::string &key){
    if (this->hist.find (key) == this->hist.end()){
        return 0.0;
    } else {
        return (this->hist.at(key));
    }
}

void CASA::loadSASAFile(const std::string fsasa){
    std::ifstream sasaFile;
    std::istream* sasainp;
    std::string line;
    std::vector<std::string> s;
    if (fsasa.length() > 0){
        sasaFile.open(fsasa.c_str(), std::ios::in);
        sasainp=&sasaFile;
        while (sasainp->good() && !(sasainp->eof())){
            getline(*sasainp, line);
            Misc::splitStr(line, " ", s, true);
            if (s.size() >= 3){
                this->targetSASA.insert(std::pair<std::string,double>(Misc::trim(s.at(0))+":"+Misc::trim(s.at(1)),atof(Misc::trim(s.at(2)).c_str())));
            }
        }
    }
}

double CASA::getTargetSASA(const std::string &key){
    if (this->targetSASA.find (key) == this->targetSASA.end()){
        return 0.0;
    } else {
        return (this->targetSASA.at(key));
    }
}

void CASA::loadParmFile(const std::string fparmfile)
{
	std::ifstream parmFile;
	std::istream* parminp;
	std::string line, key;
	std::vector<std::string> s;
	double alpha;
	if (fparmfile.length() > 0){
		parmFile.open(fparmfile.c_str(), std::ios::in);
		parminp=&parmFile;
		if (parmFile.is_open()){
			while (parminp->good() && !(parminp->eof()))
			{
				getline(*parminp, line);
				Misc::splitStr(line, " ", s, true);            
				if (s.size() ==  2){
					key = Misc::trim(s.at(0));
					alpha = atof(Misc::trim(s.at(1)).c_str());
					this->sasaCoef.insert(std::pair<std::string,double>(key,alpha));
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << fparmfile << std::endl;	
			exit(0);		
		}
	}
}

void CASA::loadRefFile(const std::string freffile)
{
	std::ifstream refFile;
	std::istream* refinp;
	std::string line, key;
	std::vector<std::string> s;
	double refsasa;
	if (freffile.length() > 0){
		refFile.open(freffile.c_str(), std::ios::in);
		refinp=&refFile;
		if (refFile.is_open()) {
			while (refinp->good() && !(refinp->eof()))
			{
				getline(*refinp, line);
				Misc::splitStr(line, " ", s, true);
				if (s.size() ==  2){
					key = Misc::trim(s.at(0));				
					refsasa = atof(Misc::trim(s.at(1)).c_str());
					this->refSASA.insert(std::pair<std::string,double>(key,refsasa));	
				}
			}
		}
		else {
			std::cerr << "Error: Unable to open file: " << freffile << std::endl;	
			exit(0);		                
		}
	}
}

