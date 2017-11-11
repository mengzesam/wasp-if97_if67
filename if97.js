var Wasp97=function(){
	/*
	reference to "The International Association for the Properties of Water and Steam:
				Revised Release on the IAPWS Industrial Formulation 1997
				for the Thermodynamic Properties of Water and Steam
				(The revision only relates to the extension of region 5 to 50 MPa)" . The following is abbrev. as if97
	*/
	/*
	p:Mpa, t:℃, T:K, h:kJ/kg, s:kJ/kg.K,v:kg/㎥
	*/
	if (arguments.length == 1 && typeof(arguments[0]) == 'number') {
		var in1 = arguments[0];
		var in_method = 'p';
	} else if (arguments.length == 2 && typeof(arguments[0]) == 'number' && typeof(arguments[1]) == 'number') {
		var in1 = arguments[0];
		var in2 = arguments[1];
		var in_method = 'pt';
	} else if (arguments.length == 2 && typeof(arguments[0]) == 'number' && typeof(arguments[1]) == 'string') {
		var in1 = arguments[0];
		var in_method = arguments[1];
		if (in_method != 'p' || in_method != 't') {
			in_method = 'p'
		}
	} else if (arguments.length == 3 && typeof(arguments[0]) == 'number' && typeof(arguments[1]) == 'number' && typeof(arguments[2]) == 'string') {
		var in1 = arguments[0];
		var in2 = arguments[1];
		var in_method = arguments[2];
		if (in_method != 'pt' && in_method != 'ph' && in_method != 'ps' && in_method != 'pv' && in_method != 'px' && in_method != 'th' && in_method != 'ts' && in_method != 'tv' && in_method != 'tx') {
			in_method = 'pt'
		}
	} else {
		var in1 = 1;
		var in2 = 200;
		var in_method = 'pt';
	}
	this.in_method=in_method;
	this.p=0;
	this.t=0;
	this.h=0;
	this.s=0;
	this.v=0;
	this.x=0;
	switch (this.in_method) {
		case 'pt':
			this.p = in1;
			this.t = in2;
			break;
		case 'ph':
			this.p = in1;
			this.h = in2;
			break;
		case 'ps':
			this.p = in1;
			this.s = in2;
			break;
		case 'pv':
			this.p = in1;
			this.v = in2;
			break;
		case 'px':
			this.p = in1;
			this.x = in2;
			break;
		case 'th':
			this.t = in1;
			this.h = in2;
			break;
		case 'ts':
			this.t = in1;
			this.s = in2;
			break;
		case 'tv':
			this.t = in1;
			this.v = in2;
			break;
		case 'tx':
			this.t = in1;
			this.x = in2;
			break;
		case 'p':
			this.p = in1;
			break;
		case 't':
			this.t = in1;
			break;
	}

	Wasp97.prototype.R=0.461526;
	Wasp97.prototype.Tc= 647.096;
	Wasp97.prototype.Pc = 22.064;
	Wasp97.prototype.rhoc = 322;
	Wasp97.prototype.T0=273.15;
	Wasp97.prototype.boundary={
		Tmin1:273.15,Tmax1:1073.15,
		Tmin2:1073.15,Tmax2:2273.15,
		Pmin1:0,Pmax1:100,
		Pmin2:0,Pmax2:50
	}
	Wasp97.prototype.table1={ //if97 table1 page6
		ni:[
			0.34805185628969E3,
			-0.11671859879975E1,
			0.10192970039326E-2,
			0.57254459862746E3,
			0.13918839778870E2
		]
	};
	Wasp97.prototype.table2={ //if97 table2 page7
		ni:[
			 0.14632971213167E0,
			-0.84548187169114E0,
			-0.37563603672040E1,
			 0.33855169168385E1,
			-0.95791963387872E0,
			 0.15772038513228E0,
			-0.16616417199501E-1,
			 0.81214629983568E-3,
			 0.28319080123804E-3,
			-0.60706301565874E-3,
			-0.18990068218419E-1,
			-0.32529748770505E-1,
			-0.21841717175414E-1,
			-0.52838357969930E-4,
			-0.47184321073267E-3,
			-0.30001780793026E-3,
			 0.47661393906987E-4,
			-0.44141845330846E-5,
			-0.72694996297594E-15,
			-0.31679644845054E-4,
			-0.28270797985312E-5,
			-0.85205128120103E-9,
			-0.22425281908000E-5,
			-0.65171222895601E-6,
			-0.14341729937924E-12,
			-0.40516996860117E-6,
			-0.12734301741641E-8,
			-0.17424871230634E-9,
			-0.68762131295531E-18,
			 0.14478307828521E-19,
			 0.26335781662795E-22,
			-0.11947622640071E-22,
			 0.18228094581404E-23,
			-0.93537087292458E-25],
		Ii:[
			0,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,3,4,4,4,5,8,8,21,23,29,30,31,32],
		Ji:[
			-2,-1,0,1,2,3,4,5,-9,-7,-1,0,1,3,-3,0,1,3,17,-4,0,6,-5,-2,10,-8,-11,-6,-29,-31,-38,-39,-40,-41]
	};
	Wasp97.prototype.table6={ //if97 table6 page10
		ni:[
			-0.23872489924521E3,
			 0.40421188637945E3,
			 0.11349746881718E3,
			-0.58457616048039E1,
			-0.15285482413140E-3,
			-0.10866707695377E-5,
			-0.13391744872602E2,
			 0.43211039183559E2,
			-0.54010067170506E2,
			 0.30535892203916E2,
			-0.65964749423638E1,
			 0.93965400878363E-2,
			 0.11573647505340E-6,
			-0.25858641282073E-4,
			-0.40644363084799E-8,
			 0.66456186191635E-7,
			 0.80670734103027E-10,
			-0.93477771213947E-12,
			 0.58265442020601E-14,
			-0.15020185953503E-16],
		Ii:[
			0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,3,3,4,5,6],
		Ji:[
			0,1,2,6,22,32,0,1,2,3,4,10,32,10,32,10,32,32,32,32]
	}
	Wasp97.prototype.table8={ //if97 table8 page11
		ni:[
			 0.17478268058307E3,
			 0.34806930892873E2,
			 0.65292584978455E1,
			 0.33039981775489E0,
			-0.19281382923196E-6,
			-0.24909197244573E-22,
			-0.26107636489332E0,
			 0.22592965981586E0,
			-0.64256463395226E-1,
			 0.78876289270526E-2,
			 0.35672110607366E-9,
			 0.17332496994895E-23,
			 0.56608900654837E-3,
			-0.32635483139717E-3,
			 0.44778286690632E-4,
			-0.51322156908507E-9,
			-0.42522657042207E-25,
			 0.26400441360689E-12,
			 0.78124600459723E-28,
			-0.30732199903668E-30],
		Ii:[
			0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,3,3,4],
		Ji:[
			0,1,2,3,11,31,0,1,2,3,12,31,0,1,2,9,31,10,32,32]
	};
	Wasp97.prototype.table10={ //if97 table10 page13
		ni:[
			-0.96927686500217E1,
			 0.10086655968018E2,
			-0.56087911283020E-2,
			 0.71452738081455E-1,
			-0.40710498223928E0,
			 0.14240819171444E1,
			-0.43839511319450E1,
			-0.28408632460772E0,
			 0.21268463753307E-1],
		Ji:[
			0,1,-5,-4,-3,-2,-1,2,3]
	};
	Wasp97.prototype.table11={ //if97 table11 page14
		ni:[
			-0.17731742473213E-2,
			-0.17834862292358E-1,
			-0.45996013696365E-1,
			-0.57581259083432E-1,
			-0.50325278727930E-1,
			-0.33032641670203E-4,
			-0.18948987516315E-3,
			-0.39392777243355E-2,
			-0.43797295650573E-1,
			-0.26674547914087E-4,
			 0.20481737692309E-7,
			 0.43870667284435E-6,
			-0.32277677238570E-4,
			-0.15033924542148E-2,
			-0.40668253562649E-1,
			-0.78847309559367E-9,
			 0.12790717852285E-7,
			 0.48225372718507E-6,
			 0.22922076337661E-5,
			-0.16714766451061E-10,
			-0.21171472321355E-2,
			-0.23895741934104E2,
			-0.59059564324270E-17,
			-0.12621808899101E-5,
			-0.38946842435739E-1,
			 0.11256211360459E-10,
			-0.82311340897998E1,
			 0.19809712802088E-7,
			 0.10406965210174E-18,
			-0.10234747095929E-12,
			-0.10018179379511E-8,
			-0.80882908646985E-10,
			 0.10693031879409E0,
			-0.33662250574171E0,
			 0.89185845355421E-24,
			 0.30629316876232E-12,
			-0.42002467698208E-5,
			-0.59056029685639E-25,
			 0.37826947613457E-5,
			-0.12768608934681E-14,
			 0.73087610595061E-28,
			 0.55414715350778E-16,
			-0.94369707241210E-6],
		Ii:[
			1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,5,6,6,6,7,7,7,8,8,9,10,10,10,16,16,18,20,20,20,21,22,23,24,24,24],
		Ji:[
			0,1,2,3,6,1,2,4,7,36,0,1,3,6,35,1,2,3,7,3,16,35,0,11,25,8,36,13,4,10,14,29,50,57,20,35,48,21,53,39,26,40,58]
	};
	Wasp97.prototype.table16={ //if97 table16 page18
		ni:[
			-0.73362260186506E-2,
			-0.88223831943146E-1,
			-0.72334555213245E-1,
			-0.40813178534455E-2,
			 0.20097803380207E-2,
			-0.53045921898642E-1,
			-0.76190409086970E-2,
			-0.63498037657313E-2,
			-0.86043093028588E-1,
			 0.75321581522770E-2,
			-0.79238375446139E-2,
			-0.22888160778447E-3,
			-0.26456501482810E-2],
		Ii:[
			1,1,1,1,2,2,2,3,3,4,4,5,5],
		Ji:[
			0,2,5,11,1,7,16,4,16,7,10,9,10]
	};
	Wasp97.prototype.table19={ //if97 table19 page22
		ni:[
			 0.90584278514723E3,
			-0.67955786399241E0,
			 0.12809002730136E-3,
			 0.26526571908428E4,
			 0.45257578905948E1]
	};
	Wasp97.prototype.table20={ //if97 table20 page22
		ni:[
			 0.10898952318288E4,
			 0.84951654495535E3,
			-0.10781748091826E3,
			 0.33153654801263E2,
			-0.74232016790248E1,
			 0.11765048724356E2,
			 0.18445749355790E1,
			-0.41792700549624E1,
			 0.62478196935812E1,
			-0.17344563108114E2,
			-0.20058176862096E3,
			 0.27196065473796E3,
			-0.45511318285818E3,
			 0.30919688604755E4,
			 0.25226640357872E6,
			-0.61707422868339E-2,
			-0.31078046629583E0,
			 0.11670873077107E2,
			 0.12812798404046E9,
			-0.98554909623276E9,
			 0.28224546973002E10,
			-0.35948971410703E10,
			 0.17227349913197E10,
			-0.13551334240775E5,
			 0.12848734664650E8,
			 0.13865724283226E1,
			 0.23598832556514E6,
			-0.13105236545054E8,
			 0.73999835474766E4,
			-0.55196697030060E6,
			 0.37154085996233E7,
			 0.19127729239660E5,
			-0.41535164835634E6,
			-0.62459855192507E2],
		Ii:[
			0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,4,4,4,5,5,5,6,6,7],
		Ji:[
			0,1,2,3,7,20,0,1,2,3,7,9,11,18,44,0,2,7,36,38,40,42,44,24,44,12,32,44,32,36,42,34,44,28]
	};
	Wasp97.prototype.table21={ //if97 table21 page23
		ni:[
			 0.14895041079516E4,
			 0.74307798314034E3,
			-0.97708318797837E2,
			 0.24742464705674E1,
			-0.63281320016026E0,
			 0.11385952129658E1,
			-0.47811863648625E0,
			 0.85208123431544E-2,
			 0.93747147377932E0,
			 0.33593118604916E1,
			 0.33809355601454E1,
			 0.16844539671904E0,
			 0.73875745236695E0,
			-0.47128737436186E0,
			 0.15020273139707E0,
			-0.21764114219750E-2,
			-0.21810755324761E-1,
			-0.10829784403677E0,
			-0.46333324635812E-1,
			 0.71280351959551E-4,
			 0.11032831789999E-3,
			 0.18955248387902E-3,
			 0.30891541160537E-2,
			 0.13555504554949E-2,
			 0.28640237477456E-6,
			-0.10779857357512E-4,
			-0.76462712454814E-4,
			 0.14052392818316E-4,
			-0.31083814331434E-4,
			-0.10302738212103E-5,
			 0.28217281635040E-6,
			 0.12704902271945E-5,
			 0.73803353468292E-7,
			-0.11030139238909E-7,
			-0.81456365207833E-13,
			-0.25180545682962E-10,
			-0.17565233969407E-17,
			 0.86934156344163E-14],
		Ii:[
			0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,4,5,5,5,6,7,7,9,9],
		Ji:[
			0,1,2,12,18,24,28,40,0,2,6,12,18,24,28,40,2,8,18,40,1,2,12,24,2,12,18,24,28,40,18,24,40,28,2,28,1,40]
	};
	Wasp97.prototype.table22={ //if97 table22 page24
		ni:[
			-0.32368398555242E13,
			 0.73263350902181E13,
			 0.35825089945447E12,
			-0.58340131851590E12,
			-0.10783068217470E11,
			 0.20825544563171E11,
			 0.61074783564516E6,
			 0.85977722535580E6,
			-0.25745723604170E5,
			 0.31081088422714E5,
			 0.12082315865936E4,
			 0.48219755109255E3,
			 0.37966001272486E1,
			-0.10842984880077E2,
			-0.45364172676660E-1,
			 0.14559115658698E-12,
			 0.11261597407230E-11,
			-0.17804982240686E-10,
			 0.12324579690832E-6,
			-0.11606921130984E-5,
			 0.27846367088554E-4,
			-0.59270038474176E-3,
			 0.12918582991878E-2],
		Ii:[
			-7,-7,-6,-6,-5,-5,-2,-2,-1,-1,0,0,1,1,2,6,6,6,6,6,6,6,6],
		Ji:[
			0,4,0,2,0,2,0,1,0,2,0,1,4,8,4,0,1,4,10,12,16,20,22]
	};
	Wasp97.prototype.table25={ //if97 table25 page26
		ni:[
			-0.39235983861984E6,
			 0.51526573827270E6,
			 0.40482443161048E5,
			-0.32193790923902E3,
			 0.96961424218694E2,
			-0.22867846371773E2,
			-0.44942914124357E6,
			-0.50118336020166E4,
			 0.35684463560015E0,
			 0.44235335848190E5,
			-0.13673388811708E5,
			 0.42163260207864E6,
			 0.22516925837475E5,
			 0.47442144865646E3,
			-0.14931130797647E3,
			-0.19781126320452E6,
			-0.23554399470760E5,
			-0.19070616302076E5,
			 0.55375669883164E5,
			 0.38293691437363E4,
			-0.60391860580567E3,
			 0.19363102620331E4,
			 0.42660643698610E4,
			-0.59780638872718E4,
			-0.70401463926862E3,
			 0.33836784107553E3,
			 0.20862786635187E2,
			 0.33834172656196E-1,
			-0.43124428414893E-4,
			 0.16653791356412E3,
			-0.13986292055898E3,
			-0.78849547999872E0,
			 0.72132411753872E-1,
			-0.59754839398283E-2,
			-0.12141358953904E-4,
			 0.23227096733871E-6,
			-0.10538463566194E2,
			 0.20718925496502E1,
			-0.72193155260427E-1,
			 0.20749887081120E-6,
			-0.18340657911379E-1,
			 0.29036272348696E-6,
			 0.21037527893619E0,
			 0.25681239729999E-3,
			-0.12799002933781E-1,
			-0.82198102652018E-5],
		Ii:[
			-1.5,-1.5,-1.5,-1.5,-1.5,-1.5,-1.25,-1.25,-1.25,-1,-1,-1,-1,-1,-1,-0.75,-0.75,-0.5,-0.5,-0.5,-0.5,-0.25,
			-0.25,-0.25,-0.25,0.25,0.25,0.25,0.25,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.75,0.75,0.75,0.75,1,1,1.25,1.25,1.5,1.5],
		Ji:[
			-24,-23,-19,-13,-11,-10,-19,-15,-6,-26,-21,-17,-16,-9,-8,-15,-14,-26,-13,-9,
			-7,-27,-25,-11,-6,1,4,8,11,0,1,5,6,10,14,16,0,4,9,17,7,18,3,15,5,18]
	};
	Wasp97.prototype.table26={ //if97 table26 page27
		ni:[
			 0.31687665083497E6,
			 0.20864175881858E2,
			-0.39859399803599E6,
			-0.21816058518877E2,
			 0.22369785194242E6,
			-0.27841703445817E4,
			 0.99207436071480E1,
			-0.75197512299157E5,
			 0.29708605951158E4,
			-0.34406878548526E1,
			 0.38815564249115E0,
			 0.17511295085750E5,
			-0.14237112854449E4,
			 0.10943803364167E1,
			 0.89971619308495E0,
			-0.33759740098958E4,
			 0.47162885818355E3,
			-0.19188241993679E1,
			 0.41078580492196E0,
			-0.33465378172097E0,
			 0.13870034777505E4,
			-0.40663326195838E3,
			 0.41727347159610E2,
			 0.21932549434532E1,
			-0.10320050009077E1,
			 0.35882943516703E0,
			 0.52511453726066E-2,
			 0.12838916450705E2,
			-0.28642437219381E1,
			 0.56912683664855E0,
			-0.99962954584931E-1,
			-0.32632037778459E-2,
			 0.23320922576723E-3,
			-0.15334809857450E0,
			 0.29072288239902E-1,
			 0.37534702741167E-3,
			 0.17296691702411E-2,
			-0.38556050844504E-3,
			-0.35017712292608E-4,
			-0.14566393631492E-4,
			 0.56420857267269E-5,
			 0.41286150074605E-7,
			-0.20684671118824E-7,
			 0.16409393674725E-8],
		Ii:[
			-6,-6,-5,-5,-4,-4,-4,-3,-3,-3,-3,-2,-2,-2,-2,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,3,3,3,4,4,5,5,5],
		Ji:[
			0,11,0,11,0,1,11,0,1,11,12,0,1,6,10,0,1,5,8,9,0,1,2,4,5,6,9,0,1,2,3,7,8,0,1,5,0,1,3,0,1,0,1,2]
	};
	Wasp97.prototype.table27={ //if97 table27 page28
		ni:[
			 0.90968501005365E3,
			 0.24045667088420E4,
			-0.59162326387130E3,
			 0.54145404128074E3,
			-0.27098308411192E3,
			 0.97976525097926E3,
			-0.46966772959435E3,
			 0.14399274604723E2,
			-0.19104204230429E2,
			 0.53299167111971E1,
			-0.21252975375934E2,
			-0.31147334413760E0,
			 0.60334840894623E0,
			-0.42764839702509E-1,
			 0.58185597255259E-2,
			-0.14597008284753E-1,
			 0.56631175631027E-2,
			-0.76155864584577E-4,
			 0.22440342919332E-3,
			-0.12561095013413E-4,
			 0.63323132660934E-6,
			-0.20541989675375E-5,
			 0.36405370390082E-7,
			-0.29759897789215E-8,
			 0.10136618529763E-7,
			 0.59925719692351E-11,
			-0.20677870105164E-10,
			-0.20874278181886E-10,
			 0.10162166825089E-9,
			-0.16429828281347E-9],
		Ii:[
			-2,-2,-1,0,0,0,0,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7,7,7,7,7],
		Ji:[
			0,1,0,0,1,2,3,0,1,3,4,0,1,2,0,1,5,0,1,4,0,1,2,0,1,0,1,3,4,5]
	};
	Wasp97.prototype.table30={ //if97 table30 page30
		ni:[
			 0.10658070028513E1,
			-0.15732845290239E2,
			 0.20944396974307E2,
			-0.76867707878716E1,
			 0.26185947787954E1,
			-0.28080781148620E1,
			 0.12053369696517E1,
			-0.84566812812502E-2,
			-0.12654315477714E1,
			-0.11524407806681E1,
			 0.88521043984318E0,
			-0.64207765181607E0,
			 0.38493460186671E0,
			-0.85214708824206E0,
			 0.48972281541877E1,
			-0.30502617256965E1,
			 0.39420536879154E-1,
			 0.12558408424308E0,
			-0.27999329698710E0,
			 0.13899799569460E1,
			-0.20189915023570E1,
			-0.82147637173963E-2,
			-0.47596035734923E0,
			 0.43984074473500E-1,
			-0.44476435428739E0,
			 0.90572070719733E0,
			 0.70522450087967E0,
			 0.10770512626332E0,
			-0.32913623258954E0,
			-0.50871062041158E0,
			-0.22175400873096E-1,
			 0.94260751665092E-1,
			 0.16436278447961E0,
			-0.13503372241348E-1,
			-0.14834345352472E-1,
			 0.57922953628084E-3,
			 0.32308904703711E-2,
			 0.80964802996215E-4,
			-0.16557679795037E-3,
			-0.44923899061815E-4],
		Ii:[
			-999,0,0,0,0,0,0,0,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,6,6,6,7,8,9,9,10,10,11],
		Ji:[
			-999,0,1,2,7,10,12,23,2,6,15,17,0,2,6,7,22,26,0,2,4,16,26,0,2,4,26,1,3,26,0,2,26,2,26,2,26,0,1,26]
	};
	Wasp97.prototype.table34={ //if97 table34 page33
		ni:[
			 0.11670521452767E4,
			-0.72421316703206E6,
			-0.17073846940092E2,
			 0.12020824702470E5,
			-0.32325550322333E7,
			 0.14915108613530E2,
			-0.48232657361591E4,
			 0.40511340542057E6,
			-0.23855557567849E0,
			 0.65017534844798E3]
	};
	Wasp97.prototype.table37={ //if97 table37 page36
		ni:[
			-0.13179983674201E2,
			 0.68540841634434E1,
			-0.24805148933466E-1,
			 0.36901534980333E0,
			-0.31161318213925E1,
			-0.32961626538917E0],
		Ji:[
			0,1,-3,-2,-1,2]
	};
	Wasp97.prototype.table38={ //if97 table38 page37
		ni:[
			 0.15736404855259E-2,
			 0.90153761673944E-3,
			-0.50270077677648E-2,
			 0.22440037409485E-5,
			-0.41163275453471E-5,
			 0.37919454822955E-7],
		Ii:[
			1,1,1,2,2,3],
		Ji:[
			1,2,3,3,9,7]
	};

	Wasp97.prototype.B23_eq5=function(t){
		var Pstar=1; //p* 1Mpa
		var Tstar=1; //T* 1K
		var theta=(t+this.T0)/Tstar;		
		var pi=this.table1.ni[0]+this.table1.ni[1]*theta+this.table1.ni[2]*Math.pow(theta,2);
		return Pstar*pi;
	}

	Wasp97.prototype.B23_eq6=function(p){
		var Pstar=1; //p* 1Mpa
		var Tstar=1; //T* 1K
		var pi=p/Pstar;
		var theta=this.table1.ni[3]+Math.sqrt((pi-this.table1.ni[4])/this.table1.ni[2]);
		return Tstar*theta-this.T0;
	}

	Wasp97.prototype.Reg1_eq7=function(p,t){//return Gibbs free energy :γ= g/(RT )
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma=0;
		for (var i = 0;i<ni.length; i++) {
			gamma+=ni[i]*Math.pow(7.1-pi,Ii[i])*Math.pow(tau-1.222,Ji[i]);
		}
		return gamma;
	}

	Wasp97.prototype.Reg1_eq7_pi=function(p,t){//return :γπ,see to if97--table4
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma_pi=0;
		for (var i = 0;i<ni.length; i++) {
			gamma_pi+=-ni[i]*Ii[i]*Math.pow(7.1-pi,Ii[i]-1)*Math.pow(tau-1.222,Ji[i]);
		}
		return gamma_pi;
	}

	Wasp97.prototype.Reg1_eq7_tau=function(p,t){//return :γτ,see to if97--table4
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma_tau=0;
		for (var i = 0;i<ni.length; i++) {
			gamma_tau+=ni[i]*Math.pow(7.1-pi,Ii[i])*Ji[i]*Math.pow(tau-1.222,Ji[i]-1);
		}
		return gamma_tau;
	}

	Wasp97.prototype.Reg1_eq7_pi2=function(p,t){//return :γππ,see to if97--table4
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma_pi2=0;
		for (var i = 0;i<ni.length; i++) {
			gamma_pi2+=ni[i]*Ii[i]*(Ii[i]-1)*Math.pow(7.1-pi,Ii[i]-2)*Math.pow(tau-1.222,Ji[i]);
		}
		return gamma_pi2;
	}

	Wasp97.prototype.Reg1_eq7_tau2=function(p,t){//return :γττ,see to if97--table4
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma_tau2=0;
		for (var i = 0;i<ni.length; i++) {
			gamma_tau2+=ni[i]*Math.pow(7.1-pi,Ii[i])*Ji[i]*(Ji[i]-1)*Math.pow(tau-1.222,Ji[i]-2);
		}
		return gamma_tau2;
	}

	Wasp97.prototype.Reg1_eq7_pitau=function(p,t){//return :γπτ,see to if97--table4
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table2.ni;
		var Ii=this.table2.Ii;
		var Ji=this.table2.Ji;
		var gamma_pitau=0;
		for (var i = 0;i<ni.length; i++) {
			gamma_pitau+=-ni[i]*Ii[i]*Math.pow(7.1-pi,Ii[i]-1)*Ji[i]*Math.pow(tau-1.222,Ji[i]-1);
		}
		return gamma_pitau;
	}

	Wasp97.prototype.Reg1_pt2h=function(p,t){//return :Specific enthalpy,see to if97--table3
		var Tstar=1386; //T* 1386K
		var tau=Tstar/(t+this.T0);
		var gamma_tau=this.Reg1_eq7_tau(p,t);
		var h=tau*gamma_tau*(t+this.T0)*this.R;
		return h;
	}

	Wasp97.prototype.Reg1_pt2s=function(p,t){//return :Specific entropy,see to if97--table3
		var Tstar=1386; //T* 1386K
		var tau=Tstar/(t+this.T0);
		var gamma=this.Reg1_eq7(p,t);
		var gamma_tau=this.Reg1_eq7_tau(p,t);
		var s=(tau*gamma_tau-gamma)*this.R;
		return s;
	}

	Wasp97.prototype.Reg1_pt2v=function(p,t){//return :Specific volume,see to if97--table3
		var Pstar=16.53; //p* 16.53Mpa
		var pi=p/Pstar;
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var v=pi*gamma_pi*(t+this.T0)*this.R/(p*1000);  //(p*1000)???
		return v;
	}

	Wasp97.prototype.Reg1_pt2u=function(p,t){//return :Specific internal energy,see to if97--table3
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var gamma_tau=this.Reg1_eq7_tau(p,t);
		var u=(tau*gamma_tau-pi*gamma_pi)*(t+this.T0)*this.R;
		return u;
	}

	Wasp97.prototype.Reg1_pt2cp=function(p,t){//return :Specific isobaric heat capacity,see to if97--table3
		var Tstar=1386; //T* 1386K
		var tau=Tstar/(t+this.T0);
		var gamma_tau2=this.Reg1_eq7_tau2(p,t);
		var cp=-tau*tau*gamma_tau2*this.R;
		return cp;
	}

	Wasp97.prototype.Reg1_pt2cv=function(p,t){//return :Specific isochoric heat capacity,see to if97--table3
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var gamma_tau2=this.Reg1_eq7_tau2(p,t);
		var gamma_pi2=this.Reg1_eq7_pi2(p,t);
		var gamma_pitau=this.Reg1_eq7_pitau(p,t);
		var cv=(-tau*tau*gamma_tau2+(gamma_pi-tau*gamma_pitau)*(gamma_pi-tau*gamma_pitau)/gamma_pi2)*this.R;
		return cv;
	}

	Wasp97.prototype.Reg1_pt2w=function(p,t){//return :Speed of sound,see to if97--table3
		var Pstar=16.53; //p* 16.53Mpa
		var Tstar=1386; //T* 1386K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var gamma_tau2=this.Reg1_eq7_tau2(p,t);
		var gamma_pi2=this.Reg1_eq7_pi2(p,t);
		var gamma_pitau=this.Reg1_eq7_pitau(p,t);
		var w=(gamma_pi*gamma_pi)/(Math.pow(gamma_pi-tau*gamma_pitau,2)/(tau*tau*gamma_pi)-gamma_pi2);
		w=w*(t+this.T0)*this.R;
		w=Math.sqrt(w);
		return w; //formula is error
	}

	Wasp97.prototype.Reg1_pt2av=function(p,t){//return :isobaric cubic expansion coefficient,see to if97--table3
		var Tstar=1386; //T* 1386K
		var tau=Tstar/(t+this.T0);
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var gamma_pitau=this.Reg1_eq7_pitau(p,t);
		var av=(1-tau*gamma_pitau/gamma_pi)/(t+this.T0);		
		return av;
	}

	Wasp97.prototype.Reg1_pt2kT=function(p,t){//return :isothermal compressibility,see to if97--table3
		var Pstar=16.53; //p* 16.53Mpa
		var pi=p/Pstar;
		var gamma_pi=this.Reg1_eq7_pi(p,t);
		var gamma_pi2=this.Reg1_eq7_pi2(p,t);
		var kT=(-pi*gamma_pi2/gamma_pi)/p;		
		return kT;
	}

	Wasp97.prototype.Reg2_eq16=function(p,t){//return Gibbs free energy ideal-gas part:γ0, g/(RT )=γ0(π,τ)+γr(π,τ)
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table10.ni;
		var Ji=this.table10.Ji;
		var gamma0=Math.log(pi);
		for (var i = 0;i<ni.length; i++) {
			gamma0+=ni[i]*Math.pow(tau,Ji[i]);
		}
		return gamma0;
	}

	Wasp97.prototype.Reg2_eq16_pi=function(p,t){//return :γ0π,see to if97--table13
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		var gamma0_pi=1/pi;
		return gamma0_pi;
	}

	Wasp97.prototype.Reg2_eq16_pi2=function(p,t){//return :γ0ππ,see to if97--table13
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		var gamma0_pi2=-1/(pi*pi);
		return gamma0_pi2;
	}

	Wasp97.prototype.Reg2_eq16_tau=function(p,t){//return :γ0τ,see to if97--table13
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var ni=this.table10.ni;
		var Ji=this.table10.Ji;
		var gamma0_tau=0;
		for (var i = 0;i<ni.length; i++) {
			gamma0_tau+=ni[i]*Ji[i]*Math.pow(tau,Ji[i]-1);
		}
		return gamma0_tau;
	}

	Wasp97.prototype.Reg2_eq16_tau2=function(p,t){//return :γ0ττ,see to if97--table13
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var ni=this.table10.ni;
		var Ji=this.table10.Ji;
		var gamma0_tau2=0;
		for (var i = 0;i<ni.length; i++) {
			gamma0_tau2+=ni[i]*Ji[i]*(Ji[i]-1)*Math.pow(tau,Ji[i]-2);
		}
		return gamma0_tau2;
	}

	Wasp97.prototype.Reg2_eq16_pitau=function(p,t){//return :γ0ττ,see to if97--table13
		var gamma0_pitau=0;
		return gamma0_pitau;
	}

	Wasp97.prototype.Reg2_eq17=function(p,t){//return Gibbs free energy residual part :γr, g/(RT )=γ0(π,τ)+γr(π,τ)
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar=0;
		for (var i = 0;i<ni.length; i++) {
			gammar+=ni[i]*Math.pow(pi,Ii[i])*Math.pow(tau-0.5,Ji[i]);
		}
		return gammar;
	}

	Wasp97.prototype.Reg2_eq17_pi=function(p,t){//return :γrπ,see to if97--table14
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar_pi=0;
		for (var i = 0;i<ni.length; i++) {
			gammar_pi+=ni[i]*Ii[i]*Math.pow(pi,Ii[i]-1)*Math.pow(tau-0.5,Ji[i]);
		}
		return gammar_pi;
	}

	Wasp97.prototype.Reg2_eq17_pi2=function(p,t){//return :γrππ,see to if97--table14
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar_pi2=0;
		for (var i = 0;i<ni.length; i++) {
			gammar_pi2+=ni[i]*Ii[i]*(Ii[i]-1)*Math.pow(pi,Ii[i]-2)*Math.pow(tau-0.5,Ji[i]);
		}
		return gammar_pi2;
	}

	Wasp97.prototype.Reg2_eq17_tau=function(p,t){//return :γrτ,see to if97--table14
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar_tau=0;
		for (var i = 0;i<ni.length; i++) {
			gammar_tau+=ni[i]*Math.pow(pi,Ii[i])*Ji[i]*Math.pow(tau-0.5,Ji[i]-1);
		}
		return gammar_tau;
	}

	Wasp97.prototype.Reg2_eq17_tau2=function(p,t){//return :γrττ,see to if97--table14
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar_tau2=0;
		for (var i = 0;i<ni.length; i++) {
			gammar_tau2+=ni[i]*Math.pow(pi,Ii[i])*Ji[i]*(Ji[i]-1)*Math.pow(tau-0.5,Ji[i]-2);
		}
		return gammar_tau2;
	}

	Wasp97.prototype.Reg2_eq17_pitau=function(p,t){//return :γrπτ,see to if97--table14
		var Pstar=1; //p* 1Mpa
		var Tstar=540; //T* 540K
		var pi=p/Pstar;
		var tau=Tstar/(t+this.T0);
		var ni=this.table11.ni;
		var Ii=this.table11.Ii;
		var Ji=this.table11.Ji;
		var gammar_pitau=0;
		for (var i = 0;i<ni.length; i++) {
			gammar_pitau+=ni[i]*Ii[i]*Math.pow(pi,Ii[i]-1)*Ji[i]*Math.pow(tau-0.5,Ji[i]-1);
		}
		return gammar_pitau;
	}

	Wasp97.prototype.Reg2_pt2h=function(p,t){//return :Specific enthalpy,see to if97--table12
		//var Pstar=1; //p* 1Mpa
		//var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0_tau=this.Reg2_eq16_tau(p,t);
		var gammar_tau=this.Reg2_eq17_tau(p,t);
		var h=tau*(gamma0_tau+gammar_tau)*(t+this.T0)*this.R;
		return h;
	}

	Wasp97.prototype.Reg2_pt2s=function(p,t){//return :Specific entropy,see to if97--table12
		//var Pstar=1; //p* 1Mpa
		//var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0=this.Reg2_eq16(p,t);
		var gammar=this.Reg2_eq17(p,t);
		var gamma0_tau=this.Reg2_eq16_tau(p,t);
		var gammar_tau=this.Reg2_eq17_tau(p,t);
		var s=(tau*(gamma0_tau+gammar_tau)-(gamma0+gammar))*this.R;
		return s;
	}

	Wasp97.prototype.Reg2_pt2v=function(p,t){//return :Specific volume,see to if97--table12
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		//var Tstar=540; //T* 540K
		//var tau=Tstar/(t+this.T0);
		var gamma0_pi=this.Reg2_eq16_pi(p,t);
		var gammar_pi=this.Reg2_eq17_pi(p,t);
		var v=(pi*(gamma0_pi+gammar_pi))*(t+this.T0)*this.R/(p*1000);  //(p*1000)???;
		return v;
	}

	Wasp97.prototype.Reg2_pt2u=function(p,t){//return :Specific internal energy,see to if97--table12
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0_pi=this.Reg2_eq16_pi(p,t);
		var gammar_pi=this.Reg2_eq17_pi(p,t);
		var gamma0_tau=this.Reg2_eq16_tau(p,t);
		var gammar_tau=this.Reg2_eq17_tau(p,t);
		var u=(tau*(gamma0_tau+gammar_tau)-pi*(gamma0_pi+gammar_pi))*(t+this.T0)*this.R;
		return u;
	}

	Wasp97.prototype.Reg2_pt2cp=function(p,t){//return :Specific isobaric heat capacity,see to if97--table12
		//var Pstar=1; //p* 1Mpa
		//var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0_tau2=this.Reg2_eq16_tau2(p,t);
		var gammar_tau2=this.Reg2_eq17_tau2(p,t);
		var cp=(-tau*tau*(gamma0_tau2+gammar_tau2))*this.R;
		return cp;
	}

	Wasp97.prototype.Reg2_pt2cv=function(p,t){//return :Specific isochoric heat capacity,see to if97--table12
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0_tau2=this.Reg2_eq16_tau2(p,t);
		var gammar_tau2=this.Reg2_eq17_tau2(p,t);
		var gammar_pi=this.Reg2_eq17_pi(p,t);
		var gammar_pi2=this.Reg2_eq17_pi2(p,t);
		var gammar_pitau=this.Reg2_eq17_pitau(p,t);
		var cv=(-tau*tau*(gamma0_tau2+gammar_tau2)-Math.pow(1+pi*gammar_pi-tau*pi*gammar_pitau,2)/(1-pi*pi*gammar_pi2))*this.R;
		return cv;
	}

	Wasp97.prototype.Reg2_pt2w=function(p,t){//return :Speed of sound,see to if97--table12
		var Pstar=1; //p* 1Mpa
		var pi=p/Pstar;
		var Tstar=540; //T* 540K
		var tau=Tstar/(t+this.T0);
		var gamma0_tau2=this.Reg2_eq16_tau2(p,t);
		var gammar_tau2=this.Reg2_eq17_tau2(p,t);
		var gammar_pi=this.Reg2_eq17_pi(p,t);
		var gammar_pi2=this.Reg2_eq17_pi2(p,t);
		var gammar_pitau=this.Reg2_eq17_pitau(p,t);
		var w=1+2*pi*gammar_pi+pi*pi*gammar_pi*gammar_pi;
		w=w/(1-pi*pi*gammar_pi2+Math.pow(1+pi*gammar_pi-tau*pi*gammar_pitau,2)/(tau*tau*(gamma0_tau2+gammar_tau2)));
		w=w*(t+this.T0)*this.R;
		return w;
	}


}




var ws=new Wasp97(100,111,'pg');


p1=3;
t1=300-273.15;
p2=80;
t2=300-273.15;
p3=3;
t3=500-273.15;
console.log(ws.Reg1_pt2v(p1,t1)+"\t"+ws.Reg1_pt2v(p2,t2)+"\t"+ws.Reg1_pt2v(p3,t3));
console.log(ws.Reg1_pt2h(p1,t1)+"\t"+ws.Reg1_pt2h(p2,t2)+"\t"+ws.Reg1_pt2h(p3,t3));
console.log(ws.Reg1_pt2u(p1,t1)+"\t"+ws.Reg1_pt2u(p2,t2)+"\t"+ws.Reg1_pt2u(p3,t3));
console.log(ws.Reg1_pt2s(p1,t1)+"\t"+ws.Reg1_pt2s(p2,t2)+"\t"+ws.Reg1_pt2s(p3,t3));
console.log(ws.Reg1_pt2cp(p1,t1)+"\t"+ws.Reg1_pt2cp(p2,t2)+"\t"+ws.Reg1_pt2cp(p3,t3));
console.log(ws.Reg1_pt2cv(p1,t1)+"\t"+ws.Reg1_pt2cv(p2,t2)+"\t"+ws.Reg1_pt2cv(p3,t3));
console.log(ws.Reg1_pt2w(p1,t1)+"\t"+ws.Reg1_pt2w(p2,t2)+"\t"+ws.Reg1_pt2w(p3,t3));
console.log(ws.Reg1_pt2av(p1,t1)+"\t"+ws.Reg1_pt2av(p2,t2)+"\t"+ws.Reg1_pt2av(p3,t3));
console.log(ws.Reg1_pt2kT(p1,t1)+"\t"+ws.Reg1_pt2kT(p2,t2)+"\t"+ws.Reg1_pt2kT(p3,t3));


p1=0.0035;
t1=300-273.15;
p2=0.0035;
t2=700-273.15;
p3=30;
t3=700-273.15;
console.log(ws.Reg2_pt2v(p1,t1)+"\t"+ws.Reg2_pt2v(p2,t2)+"\t"+ws.Reg2_pt2v(p3,t3));
console.log(ws.Reg2_pt2h(p1,t1)+"\t"+ws.Reg2_pt2h(p2,t2)+"\t"+ws.Reg2_pt2h(p3,t3));
console.log(ws.Reg2_pt2u(p1,t1)+"\t"+ws.Reg2_pt2u(p2,t2)+"\t"+ws.Reg2_pt2u(p3,t3));
console.log(ws.Reg2_pt2s(p1,t1)+"\t"+ws.Reg2_pt2s(p2,t2)+"\t"+ws.Reg2_pt2s(p3,t3));
console.log(ws.Reg2_pt2cp(p1,t1)+"\t"+ws.Reg2_pt2cp(p2,t2)+"\t"+ws.Reg2_pt2cp(p3,t3));
console.log(ws.Reg2_pt2cv(p1,t1)+"\t"+ws.Reg2_pt2cv(p2,t2)+"\t"+ws.Reg2_pt2cv(p3,t3));
console.log(ws.Reg2_pt2w(p1,t1)+"\t"+ws.Reg2_pt2w(p2,t2)+"\t"+ws.Reg2_pt2w(p3,t3));
//console.log(ws.Reg2_pt2av(p1,t1)+"\t"+ws.Reg2_pt2av(p2,t2)+"\t"+ws.Reg2_pt2av(p3,t3));
//console.log(ws.Reg2_pt2kT(p1,t1)+"\t"+ws.Reg2_pt2kT(p2,t2)+"\t"+ws.Reg2_pt2kT(p3,t3));

