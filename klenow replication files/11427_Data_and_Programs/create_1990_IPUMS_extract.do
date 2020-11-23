/* Note:  Programs will run when the do files and data files are all in the same directory */
/* Note:  Raw data (*.dat files) are downloaded directly from IPUMS website */
/* Note:  2010 data is actually pooled 2010-2012 data from the American Community Survey */ 

/******************************************************************************************/
/******************************************************************************************/
/*********************** Step 1:  Read in IPUMS Data **************************************/
/******************************************************************************************/
/******************************************************************************************/
   
set more off
cd c:\ErikMain\discrimination_growth

clear
clear matrix

infix ///
 int     year                                 1-4 ///
 byte    statefip                             5-6 ///
 float  perwt                                7-12 ///
 int     relate                              13-14 ///
 int     age                                 15-17 ///
 byte    sex                                 18 ///
 int     race                                19 ///
 int     educ                                20-21 ///
 byte    empstat                             22 ///
 byte    empstatd                            23-24 ///
 byte    labforce                            25 ///
 int     occ                                 26-28 ///
 int     occ1990                             29-31 ///
 byte    wkswork1                            32-33 ///
 byte    wkswork2                            34 ///
 byte    hrswork1                            35-36 ///
 byte    hrswork2                            37 ///
 byte    uhrswork                            38-39 ///
 long    incwage                             40-45 ///
 long    incbus                              46-51 ///
 long    incfarm                             52-57 ///
 using ipums_1990_main.dat

replace perwt=perwt/100

label var year `"Census year"'
label var statefip `"State (FIPS code)"'
label var perwt `"Person weight"'
label var relate `"Relationship to household head [general version]"'
label var age `"Age"'
label var sex `"Sex"'
label var race `"Race [general version]"'
label var educ `"Educational attainment [general version]"'
label var empstat `"Employment status [general version]"'
label var empstatd `"Employment status [detailed version]"'
label var labforce `"Labor force status"'
label var occ `"Occupation"'
label var occ1990 `"Occupation, 1990 basis"'
label var wkswork1 `"Weeks worked last year"'
label var wkswork2 `"Weeks worked last year, intervalled"'
label var hrswork1 `"Hours worked last week"'
label var hrswork2 `"Hours worked last week, intervalled"'
label var uhrswork `"Usual hours worked per week"'
label var incwage `"Wage and salary income"'
label var incbus `"Non-farm business income"'
label var incfarm `"Farm income"'

label define yearlbl 1850 `"1850"'
label define yearlbl 1860 `"1860"', add
label define yearlbl 1870 `"1870"', add
label define yearlbl 1880 `"1880"', add
label define yearlbl 1900 `"1900"', add
label define yearlbl 1910 `"1910"', add
label define yearlbl 1920 `"1920"', add
label define yearlbl 1930 `"1930"', add
label define yearlbl 1940 `"1940"', add
label define yearlbl 1950 `"1950"', add
label define yearlbl 1960 `"1960"', add
label define yearlbl 1970 `"1970"', add
label define yearlbl 1980 `"1980"', add
label define yearlbl 1990 `"1990"', add
label define yearlbl 2000 `"2000"', add
label define yearlbl 2001 `"2001"', add
label define yearlbl 2002 `"2002"', add
label define yearlbl 2003 `"2003"', add
label define yearlbl 2004 `"2004"', add
label define yearlbl 2005 `"2005"', add
label define yearlbl 2006 `"2006"', add
label define yearlbl 2007 `"2007"', add
label define yearlbl 2008 `"2008"', add
label values year yearlbl

label define statefiplbl 01 `"Alabama"'
label define statefiplbl 02 `"Alaska"', add
label define statefiplbl 04 `"Arizona"', add
label define statefiplbl 05 `"Arkansas"', add
label define statefiplbl 06 `"California"', add
label define statefiplbl 08 `"Colorado"', add
label define statefiplbl 09 `"Connecticut"', add
label define statefiplbl 10 `"Delaware"', add
label define statefiplbl 11 `"District of Columbia"', add
label define statefiplbl 12 `"Florida"', add
label define statefiplbl 13 `"Georgia"', add
label define statefiplbl 15 `"Hawaii"', add
label define statefiplbl 16 `"Idaho"', add
label define statefiplbl 17 `"Illinois"', add
label define statefiplbl 18 `"Indiana"', add
label define statefiplbl 19 `"Iowa"', add
label define statefiplbl 20 `"Kansas"', add
label define statefiplbl 21 `"Kentucky"', add
label define statefiplbl 22 `"Louisiana"', add
label define statefiplbl 23 `"Maine"', add
label define statefiplbl 24 `"Maryland"', add
label define statefiplbl 25 `"Massachusetts"', add
label define statefiplbl 26 `"Michigan"', add
label define statefiplbl 27 `"Minnesota"', add
label define statefiplbl 28 `"Mississippi"', add
label define statefiplbl 29 `"Missouri"', add
label define statefiplbl 30 `"Montana"', add
label define statefiplbl 31 `"Nebraska"', add
label define statefiplbl 32 `"Nevada"', add
label define statefiplbl 33 `"New Hampshire"', add
label define statefiplbl 34 `"New Jersey"', add
label define statefiplbl 35 `"New Mexico"', add
label define statefiplbl 36 `"New York"', add
label define statefiplbl 37 `"North Carolina"', add
label define statefiplbl 38 `"North Dakota"', add
label define statefiplbl 39 `"Ohio"', add
label define statefiplbl 40 `"Oklahoma"', add
label define statefiplbl 41 `"Oregon"', add
label define statefiplbl 42 `"Pennsylvania"', add
label define statefiplbl 44 `"Rhode Island"', add
label define statefiplbl 45 `"South Carolina"', add
label define statefiplbl 46 `"South Dakota"', add
label define statefiplbl 47 `"Tennessee"', add
label define statefiplbl 48 `"Texas"', add
label define statefiplbl 49 `"Utah"', add
label define statefiplbl 50 `"Vermont"', add
label define statefiplbl 51 `"Virginia"', add
label define statefiplbl 53 `"Washington"', add
label define statefiplbl 54 `"West Virginia"', add
label define statefiplbl 55 `"Wisconsin"', add
label define statefiplbl 56 `"Wyoming"', add
label define statefiplbl 61 `"Maine-New Hampshire-Vermont"', add
label define statefiplbl 62 `"Massachusetts-Rhode Island"', add
label define statefiplbl 63 `"Minnesota-Iowa-Missouri-Kansas-Nebraska-S.Dakota-N.Dakota"', add
label define statefiplbl 64 `"Maryland-Delaware"', add
label define statefiplbl 65 `"Montana-Idaho-Wyoming"', add
label define statefiplbl 66 `"Utah-Nevada"', add
label define statefiplbl 67 `"Arizona-New Mexico"', add
label define statefiplbl 68 `"Alaska-Hawaii"', add
label define statefiplbl 72 `"Puerto Rico"', add
label define statefiplbl 97 `"Military/Mil. Reservation"', add
label define statefiplbl 99 `"State not identified"', add
label values statefip statefiplbl

label define relatelbl 01 `"Head/Householder"', add
label define relatelbl 02 `"Spouse"', add
label define relatelbl 03 `"Child"', add
label define relatelbl 04 `"Child-in-law"', add
label define relatelbl 05 `"Parent"', add
label define relatelbl 06 `"Parent-in-Law"', add
label define relatelbl 07 `"Sibling"', add
label define relatelbl 08 `"Sibling-in-Law"', add
label define relatelbl 09 `"Grandchild"', add
label define relatelbl 10 `"Other relatives"', add
label define relatelbl 11 `"Partner, friend, visitor"', add
label define relatelbl 12 `"Other non-relatives"', add
label define relatelbl 13 `"Institutional inmates"', add
label values relate relatelbl

label define agelbl 000 `"Less than 1 year old"'
label define agelbl 001 `"1"', add
label define agelbl 002 `"2"', add
label define agelbl 003 `"3"', add
label define agelbl 004 `"4"', add
label define agelbl 005 `"5"', add
label define agelbl 006 `"6"', add
label define agelbl 007 `"7"', add
label define agelbl 008 `"8"', add
label define agelbl 009 `"9"', add
label define agelbl 010 `"10"', add
label define agelbl 011 `"11"', add
label define agelbl 012 `"12"', add
label define agelbl 013 `"13"', add
label define agelbl 014 `"14"', add
label define agelbl 015 `"15"', add
label define agelbl 016 `"16"', add
label define agelbl 017 `"17"', add
label define agelbl 018 `"18"', add
label define agelbl 019 `"19"', add
label define agelbl 020 `"20"', add
label define agelbl 021 `"21"', add
label define agelbl 022 `"22"', add
label define agelbl 023 `"23"', add
label define agelbl 024 `"24"', add
label define agelbl 025 `"25"', add
label define agelbl 026 `"26"', add
label define agelbl 027 `"27"', add
label define agelbl 028 `"28"', add
label define agelbl 029 `"29"', add
label define agelbl 030 `"30"', add
label define agelbl 031 `"31"', add
label define agelbl 032 `"32"', add
label define agelbl 033 `"33"', add
label define agelbl 034 `"34"', add
label define agelbl 035 `"35"', add
label define agelbl 036 `"36"', add
label define agelbl 037 `"37"', add
label define agelbl 038 `"38"', add
label define agelbl 039 `"39"', add
label define agelbl 040 `"40"', add
label define agelbl 041 `"41"', add
label define agelbl 042 `"42"', add
label define agelbl 043 `"43"', add
label define agelbl 044 `"44"', add
label define agelbl 045 `"45"', add
label define agelbl 046 `"46"', add
label define agelbl 047 `"47"', add
label define agelbl 048 `"48"', add
label define agelbl 049 `"49"', add
label define agelbl 050 `"50"', add
label define agelbl 051 `"51"', add
label define agelbl 052 `"52"', add
label define agelbl 053 `"53"', add
label define agelbl 054 `"54"', add
label define agelbl 055 `"55"', add
label define agelbl 056 `"56"', add
label define agelbl 057 `"57"', add
label define agelbl 058 `"58"', add
label define agelbl 059 `"59"', add
label define agelbl 060 `"60"', add
label define agelbl 061 `"61"', add
label define agelbl 062 `"62"', add
label define agelbl 063 `"63"', add
label define agelbl 064 `"64"', add
label define agelbl 065 `"65"', add
label define agelbl 066 `"66"', add
label define agelbl 067 `"67"', add
label define agelbl 068 `"68"', add
label define agelbl 069 `"69"', add
label define agelbl 070 `"70"', add
label define agelbl 071 `"71"', add
label define agelbl 072 `"72"', add
label define agelbl 073 `"73"', add
label define agelbl 074 `"74"', add
label define agelbl 075 `"75"', add
label define agelbl 076 `"76"', add
label define agelbl 077 `"77"', add
label define agelbl 078 `"78"', add
label define agelbl 079 `"79"', add
label define agelbl 080 `"80"', add
label define agelbl 081 `"81"', add
label define agelbl 082 `"82"', add
label define agelbl 083 `"83"', add
label define agelbl 084 `"84"', add
label define agelbl 085 `"85"', add
label define agelbl 086 `"86"', add
label define agelbl 087 `"87"', add
label define agelbl 088 `"88"', add
label define agelbl 089 `"89"', add
label define agelbl 090 `"90 (90+ in 1980 and 1990)"', add
label define agelbl 091 `"91"', add
label define agelbl 092 `"92"', add
label define agelbl 093 `"93"', add
label define agelbl 094 `"94"', add
label define agelbl 095 `"95"', add
label define agelbl 096 `"96"', add
label define agelbl 097 `"97"', add
label define agelbl 098 `"98"', add
label define agelbl 099 `"99"', add
label define agelbl 100 `"100 (100+ in 1970)"', add
label define agelbl 101 `"101"', add
label define agelbl 102 `"102"', add
label define agelbl 103 `"103"', add
label define agelbl 104 `"104"', add
label define agelbl 105 `"105"', add
label define agelbl 106 `"106"', add
label define agelbl 107 `"107"', add
label define agelbl 108 `"108"', add
label define agelbl 109 `"109"', add
label define agelbl 110 `"110"', add
label define agelbl 111 `"111"', add
label define agelbl 112 `"112 (112+ in the 1980 internal data)"', add
label define agelbl 113 `"113"', add
label define agelbl 114 `"114"', add
label define agelbl 115 `"115 (115+ in the 1990 internal data)"', add
label define agelbl 116 `"116"', add
label define agelbl 117 `"117"', add
label define agelbl 118 `"118"', add
label define agelbl 119 `"119"', add
label define agelbl 120 `"120"', add
label define agelbl 121 `"121"', add
label define agelbl 122 `"122"', add
label define agelbl 123 `"123"', add
label define agelbl 124 `"124"', add
label define agelbl 125 `"125"', add
label define agelbl 126 `"126"', add
label define agelbl 129 `"129"', add
label define agelbl 130 `"130"', add
label define agelbl 135 `"135"', add
label values age agelbl

label define sexlbl 1 `"Male"'
label define sexlbl 2 `"Female"', add
label values sex sexlbl

label define racelbl 1 `"White"'
label define racelbl 2 `"Black/Negro"', add
label define racelbl 3 `"American Indian or Alaska Native"', add
label define racelbl 4 `"Chinese"', add
label define racelbl 5 `"Japanese"', add
label define racelbl 6 `"Other Asian or Pacific Islander"', add
label define racelbl 7 `"Other race, nec"', add
label define racelbl 8 `"Two major races"', add
label define racelbl 9 `"Three or more major races"', add
label values race racelbl

label define educlbl 00 `"N/A or no schooling"'
label define educlbl 01 `"Nursery school to grade 4"', add
label define educlbl 02 `"Grade 5, 6, 7, or 8"', add
label define educlbl 03 `"Grade 9"', add
label define educlbl 04 `"Grade 10"', add
label define educlbl 05 `"Grade 11"', add
label define educlbl 06 `"Grade 12"', add
label define educlbl 07 `"1 year of college"', add
label define educlbl 08 `"2 years of college"', add
label define educlbl 09 `"3 years of college"', add
label define educlbl 10 `"4 years of college"', add
label define educlbl 11 `"5+ years of college"', add
label values educ educlbl

label define empstatlbl 0 `"N/A"'
label define empstatlbl 1 `"Employed"', add
label define empstatlbl 2 `"Unemployed"', add
label define empstatlbl 3 `"Not in labor force"', add
label values empstat empstatlbl

label define empstatdlbl 00 `"N/A"'
label define empstatdlbl 10 `"At work"', add
label define empstatdlbl 11 `"At work, public emerg"', add
label define empstatdlbl 12 `"Has job, not working"', add
label define empstatdlbl 13 `"Armed forces"', add
label define empstatdlbl 14 `"Armed forces--at work"', add
label define empstatdlbl 15 `"Armed forces--not at work but with job"', add
label define empstatdlbl 20 `"Unemployed"', add
label define empstatdlbl 21 `"Unemp, exper worker"', add
label define empstatdlbl 22 `"Unemp, new worker"', add
label define empstatdlbl 30 `"Not in Labor Force"', add
label define empstatdlbl 31 `"NILF, housework"', add
label define empstatdlbl 32 `"NILF, unable to work"', add
label define empstatdlbl 33 `"NILF, school"', add
label define empstatdlbl 34 `"NILF, other"', add
label values empstatd empstatdlbl

label define labforcelbl 0 `"N/A"'
label define labforcelbl 1 `"No, not in the labor force"', add
label define labforcelbl 2 `"Yes, in the labor force"', add
label values labforce labforcelbl

label define occ1990lbl 003 `"Legislators"', add
label define occ1990lbl 004 `"Chief executives and public administrators"', add
label define occ1990lbl 007 `"Financial managers"', add
label define occ1990lbl 008 `"Human resources and labor relations managers"', add
label define occ1990lbl 013 `"Managers and specialists in marketing, advertising, and public relations"', add
label define occ1990lbl 014 `"Managers in education and related fields"', add
label define occ1990lbl 015 `"Managers of medicine and health occupations"', add
label define occ1990lbl 016 `"Postmasters and mail superintendents"', add
label define occ1990lbl 017 `"Managers of food-serving and lodging establishments"', add
label define occ1990lbl 018 `"Managers of properties and real estate"', add
label define occ1990lbl 019 `"Funeral directors"', add
label define occ1990lbl 021 `"Managers of service organizations, n.e.c."', add
label define occ1990lbl 022 `"Managers and administrators, n.e.c."', add
label define occ1990lbl 023 `"Accountants and auditors"', add
label define occ1990lbl 024 `"Insurance underwriters"', add
label define occ1990lbl 025 `"Other financial specialists"', add
label define occ1990lbl 026 `"Management analysts"', add
label define occ1990lbl 027 `"Personnel, HR, training, and labor relations specialists"', add
label define occ1990lbl 028 `"Purchasing agents and buyers, of farm products"', add
label define occ1990lbl 029 `"Buyers, wholesale and retail trade"', add
label define occ1990lbl 033 `"Purchasing managers, agents and buyers, n.e.c."', add
label define occ1990lbl 034 `"Business and promotion agents"', add
label define occ1990lbl 035 `"Construction inspectors"', add
label define occ1990lbl 036 `"Inspectors and compliance officers, outside construction"', add
label define occ1990lbl 037 `"Management support occupations"', add
label define occ1990lbl 043 `"Architects"', add
label define occ1990lbl 044 `"Aerospace engineer"', add
label define occ1990lbl 045 `"Metallurgical and materials engineers, variously phrased"', add
label define occ1990lbl 047 `"Petroleum, mining, and geological engineers"', add
label define occ1990lbl 048 `"Chemical engineers"', add
label define occ1990lbl 053 `"Civil engineers"', add
label define occ1990lbl 055 `"Electrical engineer"', add
label define occ1990lbl 056 `"Industrial engineers"', add
label define occ1990lbl 057 `"Mechanical engineers"', add
label define occ1990lbl 059 `"Not-elsewhere-classified engineers"', add
label define occ1990lbl 064 `"Computer systems analysts and computer scientists"', add
label define occ1990lbl 065 `"Operations and systems researchers and analysts"', add
label define occ1990lbl 066 `"Actuaries"', add
label define occ1990lbl 067 `"Statisticians"', add
label define occ1990lbl 068 `"Mathematicians and mathematical scientists"', add
label define occ1990lbl 069 `"Physicists and astronomers"', add
label define occ1990lbl 073 `"Chemists"', add
label define occ1990lbl 074 `"Atmospheric and space scientists"', add
label define occ1990lbl 075 `"Geologists"', add
label define occ1990lbl 076 `"Physical scientists, n.e.c."', add
label define occ1990lbl 077 `"Agricultural and food scientists"', add
label define occ1990lbl 078 `"Biological scientists"', add
label define occ1990lbl 079 `"Foresters and conservation scientists"', add
label define occ1990lbl 083 `"Medical scientists"', add
label define occ1990lbl 084 `"Physicians"', add
label define occ1990lbl 085 `"Dentists"', add
label define occ1990lbl 086 `"Veterinarians"', add
label define occ1990lbl 087 `"Optometrists"', add
label define occ1990lbl 088 `"Podiatrists"', add
label define occ1990lbl 089 `"Other health and therapy"', add
label define occ1990lbl 095 `"Registered nurses"', add
label define occ1990lbl 096 `"Pharmacists"', add
label define occ1990lbl 097 `"Dietitians and nutritionists"', add
label define occ1990lbl 098 `"Respiratory therapists"', add
label define occ1990lbl 099 `"Occupational therapists"', add
label define occ1990lbl 103 `"Physical therapists"', add
label define occ1990lbl 104 `"Speech therapists"', add
label define occ1990lbl 105 `"Therapists, n.e.c."', add
label define occ1990lbl 106 `"Physicians' assistants"', add
label define occ1990lbl 113 `"Earth, environmental, and marine science instructors"', add
label define occ1990lbl 114 `"Biological science instructors"', add
label define occ1990lbl 115 `"Chemistry instructors"', add
label define occ1990lbl 116 `"Physics instructors"', add
label define occ1990lbl 118 `"Psychology instructors"', add
label define occ1990lbl 119 `"Economics instructors"', add
label define occ1990lbl 123 `"History instructors"', add
label define occ1990lbl 125 `"Sociology instructors"', add
label define occ1990lbl 127 `"Engineering instructors"', add
label define occ1990lbl 128 `"Math instructors"', add
label define occ1990lbl 139 `"Education instructors"', add
label define occ1990lbl 145 `"Law instructors"', add
label define occ1990lbl 147 `"Theology instructors"', add
label define occ1990lbl 149 `"Home economics instructors"', add
label define occ1990lbl 150 `"Humanities profs/instructors, college, nec"', add
label define occ1990lbl 154 `"Subject instructors (HS/college)"', add
label define occ1990lbl 155 `"Kindergarten and earlier school teachers"', add
label define occ1990lbl 156 `"Primary school teachers"', add
label define occ1990lbl 157 `"Secondary school teachers"', add
label define occ1990lbl 158 `"Special education teachers"', add
label define occ1990lbl 159 `"Teachers , n.e.c."', add
label define occ1990lbl 163 `"Vocational and educational counselors"', add
label define occ1990lbl 164 `"Librarians"', add
label define occ1990lbl 165 `"Archivists and curators"', add
label define occ1990lbl 166 `"Economists, market researchers, and survey researchers"', add
label define occ1990lbl 167 `"Psychologists"', add
label define occ1990lbl 168 `"Sociologists"', add
label define occ1990lbl 169 `"Social scientists, n.e.c."', add
label define occ1990lbl 173 `"Urban and regional planners"', add
label define occ1990lbl 174 `"Social workers"', add
label define occ1990lbl 175 `"Recreation workers"', add
label define occ1990lbl 176 `"Clergy and religious workers"', add
label define occ1990lbl 178 `"Lawyers"', add
label define occ1990lbl 179 `"Judges"', add
label define occ1990lbl 183 `"Writers and authors"', add
label define occ1990lbl 184 `"Technical writers"', add
label define occ1990lbl 185 `"Designers"', add
label define occ1990lbl 186 `"Musician or composer"', add
label define occ1990lbl 187 `"Actors, directors, producers"', add
label define occ1990lbl 188 `"Art makers: painters, sculptors, craft-artists, and print-makers"', add
label define occ1990lbl 189 `"Photographers"', add
label define occ1990lbl 193 `"Dancers"', add
label define occ1990lbl 194 `"Art/entertainment performers and related"', add
label define occ1990lbl 195 `"Editors and reporters"', add
label define occ1990lbl 198 `"Announcers"', add
label define occ1990lbl 199 `"Athletes, sports instructors, and officials"', add
label define occ1990lbl 200 `"Professionals, n.e.c."', add
label define occ1990lbl 203 `"Clinical laboratory technologies and technicians"', add
label define occ1990lbl 204 `"Dental hygenists"', add
label define occ1990lbl 205 `"Health record tech specialists"', add
label define occ1990lbl 206 `"Radiologic tech specialists"', add
label define occ1990lbl 207 `"Licensed practical nurses"', add
label define occ1990lbl 208 `"Health technologists and technicians, n.e.c."', add
label define occ1990lbl 213 `"Electrical and electronic (engineering) technicians"', add
label define occ1990lbl 214 `"Engineering technicians, n.e.c."', add
label define occ1990lbl 215 `"Mechanical engineering technicians"', add
label define occ1990lbl 217 `"Drafters"', add
label define occ1990lbl 218 `"Surveyors, cartographers, mapping scientists and technicians"', add
label define occ1990lbl 223 `"Biological technicians"', add
label define occ1990lbl 224 `"Chemical technicians"', add
label define occ1990lbl 225 `"Other science technicians"', add
label define occ1990lbl 226 `"Airplane pilots and navigators"', add
label define occ1990lbl 227 `"Air traffic controllers"', add
label define occ1990lbl 228 `"Broadcast equipment operators"', add
label define occ1990lbl 229 `"Computer software developers"', add
label define occ1990lbl 233 `"Programmers of numerically controlled machine tools"', add
label define occ1990lbl 234 `"Legal assistants, paralegals, legal support, etc"', add
label define occ1990lbl 235 `"Technicians, n.e.c."', add
label define occ1990lbl 243 `"Supervisors and proprietors of sales jobs"', add
label define occ1990lbl 253 `"Insurance sales occupations"', add
label define occ1990lbl 254 `"Real estate sales occupations"', add
label define occ1990lbl 255 `"Financial services sales occupations"', add
label define occ1990lbl 256 `"Advertising and related sales jobs"', add
label define occ1990lbl 258 `"Sales engineers"', add
label define occ1990lbl 274 `"Salespersons, n.e.c."', add
label define occ1990lbl 275 `"Retail sales clerks"', add
label define occ1990lbl 276 `"Cashiers"', add
label define occ1990lbl 277 `"Door-to-door sales, street sales, and news vendors"', add
label define occ1990lbl 283 `"Sales demonstrators / promoters / models"', add
label define occ1990lbl 303 `"Office supervisors"', add
label define occ1990lbl 308 `"Computer and peripheral equipment operators"', add
label define occ1990lbl 313 `"Secretaries"', add
label define occ1990lbl 314 `"Stenographers"', add
label define occ1990lbl 315 `"Typists"', add
label define occ1990lbl 316 `"Interviewers, enumerators, and surveyors"', add
label define occ1990lbl 317 `"Hotel clerks"', add
label define occ1990lbl 318 `"Transportation ticket and reservation agents"', add
label define occ1990lbl 319 `"Receptionists"', add
label define occ1990lbl 323 `"Information clerks, nec"', add
label define occ1990lbl 326 `"Correspondence and order clerks"', add
label define occ1990lbl 328 `"Human resources clerks, except payroll and timekeeping"', add
label define occ1990lbl 329 `"Library assistants"', add
label define occ1990lbl 335 `"File clerks"', add
label define occ1990lbl 336 `"Records clerks"', add
label define occ1990lbl 337 `"Bookkeepers and accounting and auditing clerks"', add
label define occ1990lbl 338 `"Payroll and timekeeping clerks"', add
label define occ1990lbl 343 `"Cost and rate clerks (financial records processing)"', add
label define occ1990lbl 344 `"Billing clerks and related financial records processing"', add
label define occ1990lbl 345 `"Duplication machine operators / office machine operators"', add
label define occ1990lbl 346 `"Mail and paper handlers"', add
label define occ1990lbl 347 `"Office machine operators, n.e.c."', add
label define occ1990lbl 348 `"Telephone operators"', add
label define occ1990lbl 349 `"Other telecom operators"', add
label define occ1990lbl 354 `"Postal clerks, excluding mail carriers"', add
label define occ1990lbl 355 `"Mail carriers for postal service"', add
label define occ1990lbl 356 `"Mail clerks, outside of post office"', add
label define occ1990lbl 357 `"Messengers"', add
label define occ1990lbl 359 `"Dispatchers"', add
label define occ1990lbl 361 `"Inspectors, n.e.c."', add
label define occ1990lbl 364 `"Shipping and receiving clerks"', add
label define occ1990lbl 365 `"Stock and inventory clerks"', add
label define occ1990lbl 366 `"Meter readers"', add
label define occ1990lbl 368 `"Weighers, measurers, and checkers"', add
label define occ1990lbl 373 `"Material recording, scheduling, production, planning, and expediting clerks"', add
label define occ1990lbl 375 `"Insurance adjusters, examiners, and investigators"', add
label define occ1990lbl 376 `"Customer service reps, investigators and adjusters, except insurance"', add
label define occ1990lbl 377 `"Eligibility clerks for government programs; social welfare"', add
label define occ1990lbl 378 `"Bill and account collectors"', add
label define occ1990lbl 379 `"General office clerks"', add
label define occ1990lbl 383 `"Bank tellers"', add
label define occ1990lbl 384 `"Proofreaders"', add
label define occ1990lbl 385 `"Data entry keyers"', add
label define occ1990lbl 386 `"Statistical clerks"', add
label define occ1990lbl 387 `"Teacher's aides"', add
label define occ1990lbl 389 `"Administrative support jobs, n.e.c."', add
label define occ1990lbl 405 `"Housekeepers, maids, butlers, stewards, and lodging quarters cleaners"', add
label define occ1990lbl 407 `"Private household cleaners and servants"', add
label define occ1990lbl 415 `"Supervisors of guards"', add
label define occ1990lbl 417 `"Fire fighting, prevention, and inspection"', add
label define occ1990lbl 418 `"Police, detectives, and private investigators"', add
label define occ1990lbl 423 `"Other law enforcement: sheriffs, bailiffs, correctional institution officers"', add
label define occ1990lbl 425 `"Crossing guards and bridge tenders"', add
label define occ1990lbl 426 `"Guards, watchmen, doorkeepers"', add
label define occ1990lbl 427 `"Protective services, n.e.c."', add
label define occ1990lbl 434 `"Bartenders"', add
label define occ1990lbl 435 `"Waiter/waitress"', add
label define occ1990lbl 436 `"Cooks, variously defined"', add
label define occ1990lbl 438 `"Food counter and fountain workers"', add
label define occ1990lbl 439 `"Kitchen workers"', add
label define occ1990lbl 443 `"Waiter's assistant"', add
label define occ1990lbl 444 `"Misc food prep workers"', add
label define occ1990lbl 445 `"Dental assistants"', add
label define occ1990lbl 446 `"Health aides, except nursing"', add
label define occ1990lbl 447 `"Nursing aides, orderlies, and attendants"', add
label define occ1990lbl 448 `"Supervisors of cleaning and building service"', add
label define occ1990lbl 453 `"Janitors"', add
label define occ1990lbl 454 `"Elevator operators"', add
label define occ1990lbl 455 `"Pest control occupations"', add
label define occ1990lbl 456 `"Supervisors of personal service jobs, n.e.c."', add
label define occ1990lbl 457 `"Barbers"', add
label define occ1990lbl 458 `"Hairdressers and cosmetologists"', add
label define occ1990lbl 459 `"Recreation facility attendants"', add
label define occ1990lbl 461 `"Guides"', add
label define occ1990lbl 462 `"Ushers"', add
label define occ1990lbl 463 `"Public transportation attendants and inspectors"', add
label define occ1990lbl 464 `"Baggage porters"', add
label define occ1990lbl 465 `"Welfare service aides"', add
label define occ1990lbl 468 `"Child care workers"', add
label define occ1990lbl 469 `"Personal service occupations, nec"', add
label define occ1990lbl 473 `"Farmers (owners and tenants)"', add
label define occ1990lbl 474 `"Horticultural specialty farmers"', add
label define occ1990lbl 475 `"Farm managers, except for horticultural farms"', add
label define occ1990lbl 476 `"Managers of horticultural specialty farms"', add
label define occ1990lbl 479 `"Farm workers"', add
label define occ1990lbl 483 `"Marine life cultivation workers"', add
label define occ1990lbl 484 `"Nursery farming workers"', add
label define occ1990lbl 485 `"Supervisors of agricultural occupations"', add
label define occ1990lbl 486 `"Gardeners and groundskeepers"', add
label define occ1990lbl 487 `"Animal caretakers except on farms"', add
label define occ1990lbl 488 `"Graders and sorters of agricultural products"', add
label define occ1990lbl 489 `"Inspectors of agricultural products"', add
label define occ1990lbl 496 `"Timber, logging, and forestry workers"', add
label define occ1990lbl 498 `"Fishers, hunters, and kindred"', add
label define occ1990lbl 503 `"Supervisors of mechanics and repairers"', add
label define occ1990lbl 505 `"Automobile mechanics"', add
label define occ1990lbl 507 `"Bus, truck, and stationary engine mechanics"', add
label define occ1990lbl 508 `"Aircraft mechanics"', add
label define occ1990lbl 509 `"Small engine repairers"', add
label define occ1990lbl 514 `"Auto body repairers"', add
label define occ1990lbl 516 `"Heavy equipment and farm equipment mechanics"', add
label define occ1990lbl 518 `"Industrial machinery repairers"', add
label define occ1990lbl 519 `"Machinery maintenance occupations"', add
label define occ1990lbl 523 `"Repairers of industrial electrical equipment"', add
label define occ1990lbl 525 `"Repairers of data processing equipment"', add
label define occ1990lbl 526 `"Repairers of household appliances and power tools"', add
label define occ1990lbl 527 `"Telecom and line installers and repairers"', add
label define occ1990lbl 533 `"Repairers of electrical equipment, n.e.c."', add
label define occ1990lbl 534 `"Heating, air conditioning, and refigeration mechanics"', add
label define occ1990lbl 535 `"Precision makers, repairers, and smiths"', add
label define occ1990lbl 536 `"Locksmiths and safe repairers"', add
label define occ1990lbl 538 `"Office machine repairers and mechanics"', add
label define occ1990lbl 539 `"Repairers of mechanical controls and valves"', add
label define occ1990lbl 543 `"Elevator installers and repairers"', add
label define occ1990lbl 544 `"Millwrights"', add
label define occ1990lbl 549 `"Mechanics and repairers, n.e.c."', add
label define occ1990lbl 558 `"Supervisors of construction work"', add
label define occ1990lbl 563 `"Masons, tilers, and carpet installers"', add
label define occ1990lbl 567 `"Carpenters"', add
label define occ1990lbl 573 `"Drywall installers"', add
label define occ1990lbl 575 `"Electricians"', add
label define occ1990lbl 577 `"Electric power installers and repairers"', add
label define occ1990lbl 579 `"Painters, construction and maintenance"', add
label define occ1990lbl 583 `"Paperhangers"', add
label define occ1990lbl 584 `"Plasterers"', add
label define occ1990lbl 585 `"Plumbers, pipe fitters, and steamfitters"', add
label define occ1990lbl 588 `"Concrete and cement workers"', add
label define occ1990lbl 589 `"Glaziers"', add
label define occ1990lbl 593 `"Insulation workers"', add
label define occ1990lbl 594 `"Paving, surfacing, and tamping equipment operators"', add
label define occ1990lbl 595 `"Roofers and slaters"', add
label define occ1990lbl 596 `"Sheet metal duct installers"', add
label define occ1990lbl 597 `"Structural metal workers"', add
label define occ1990lbl 598 `"Drillers of earth"', add
label define occ1990lbl 599 `"Construction trades, n.e.c."', add
label define occ1990lbl 614 `"Drillers of oil wells"', add
label define occ1990lbl 615 `"Explosives workers"', add
label define occ1990lbl 616 `"Miners"', add
label define occ1990lbl 617 `"Other mining occupations"', add
label define occ1990lbl 628 `"Production supervisors or foremen"', add
label define occ1990lbl 634 `"Tool and die makers and die setters"', add
label define occ1990lbl 637 `"Machinists"', add
label define occ1990lbl 643 `"Boilermakers"', add
label define occ1990lbl 644 `"Precision grinders and filers"', add
label define occ1990lbl 645 `"Patternmakers and model makers"', add
label define occ1990lbl 646 `"Lay-out workers"', add
label define occ1990lbl 649 `"Engravers"', add
label define occ1990lbl 653 `"Tinsmiths, coppersmiths, and sheet metal workers"', add
label define occ1990lbl 657 `"Cabinetmakers and bench carpenters"', add
label define occ1990lbl 658 `"Furniture and wood finishers"', add
label define occ1990lbl 659 `"Other precision woodworkers"', add
label define occ1990lbl 666 `"Dressmakers and seamstresses"', add
label define occ1990lbl 667 `"Tailors"', add
label define occ1990lbl 668 `"Upholsterers"', add
label define occ1990lbl 669 `"Shoe repairers"', add
label define occ1990lbl 674 `"Other precision apparel and fabric workers"', add
label define occ1990lbl 675 `"Hand molders and shapers, except jewelers"', add
label define occ1990lbl 677 `"Optical goods workers"', add
label define occ1990lbl 678 `"Dental laboratory and medical appliance technicians"', add
label define occ1990lbl 679 `"Bookbinders"', add
label define occ1990lbl 684 `"Other precision and craft workers"', add
label define occ1990lbl 686 `"Butchers and meat cutters"', add
label define occ1990lbl 687 `"Bakers"', add
label define occ1990lbl 688 `"Batch food makers"', add
label define occ1990lbl 693 `"Adjusters and calibrators"', add
label define occ1990lbl 694 `"Water and sewage treatment plant operators"', add
label define occ1990lbl 695 `"Power plant operators"', add
label define occ1990lbl 696 `"Plant and system operators, stationary engineers"', add
label define occ1990lbl 699 `"Other plant and system operators"', add
label define occ1990lbl 703 `"Lathe, milling, and turning machine operatives"', add
label define occ1990lbl 706 `"Punching and stamping press operatives"', add
label define occ1990lbl 707 `"Rollers, roll hands, and finishers of metal"', add
label define occ1990lbl 708 `"Drilling and boring machine operators"', add
label define occ1990lbl 709 `"Grinding, abrading, buffing, and polishing workers"', add
label define occ1990lbl 713 `"Forge and hammer operators"', add
label define occ1990lbl 717 `"Fabricating machine operators, n.e.c."', add
label define occ1990lbl 719 `"Molders, and casting machine operators"', add
label define occ1990lbl 723 `"Metal platers"', add
label define occ1990lbl 724 `"Heat treating equipment operators"', add
label define occ1990lbl 726 `"Wood lathe, routing, and planing machine operators"', add
label define occ1990lbl 727 `"Sawing machine operators and sawyers"', add
label define occ1990lbl 728 `"Shaping and joining machine operator (woodworking)"', add
label define occ1990lbl 729 `"Nail and tacking machine operators  (woodworking)"', add
label define occ1990lbl 733 `"Other woodworking machine operators"', add
label define occ1990lbl 734 `"Printing machine operators, n.e.c."', add
label define occ1990lbl 735 `"Photoengravers and lithographers"', add
label define occ1990lbl 736 `"Typesetters and compositors"', add
label define occ1990lbl 738 `"Winding and twisting textile/apparel operatives"', add
label define occ1990lbl 739 `"Knitters, loopers, and toppers textile operatives"', add
label define occ1990lbl 743 `"Textile cutting machine operators"', add
label define occ1990lbl 744 `"Textile sewing machine operators"', add
label define occ1990lbl 745 `"Shoemaking machine operators"', add
label define occ1990lbl 747 `"Pressing machine operators (clothing)"', add
label define occ1990lbl 748 `"Laundry workers"', add
label define occ1990lbl 749 `"Misc textile machine operators"', add
label define occ1990lbl 753 `"Cementing and gluing maching operators"', add
label define occ1990lbl 754 `"Packers, fillers, and wrappers"', add
label define occ1990lbl 755 `"Extruding and forming machine operators"', add
label define occ1990lbl 756 `"Mixing and blending machine operatives"', add
label define occ1990lbl 757 `"Separating, filtering, and clarifying machine operators"', add
label define occ1990lbl 759 `"Painting machine operators"', add
label define occ1990lbl 763 `"Roasting and baking machine operators (food)"', add
label define occ1990lbl 764 `"Washing, cleaning, and pickling machine operators"', add
label define occ1990lbl 765 `"Paper folding machine operators"', add
label define occ1990lbl 766 `"Furnace, kiln, and oven operators, apart from food"', add
label define occ1990lbl 768 `"Crushing and grinding machine operators"', add
label define occ1990lbl 769 `"Slicing and cutting machine operators"', add
label define occ1990lbl 773 `"Motion picture projectionists"', add
label define occ1990lbl 774 `"Photographic process workers"', add
label define occ1990lbl 779 `"Machine operators, n.e.c."', add
label define occ1990lbl 783 `"Welders and metal cutters"', add
label define occ1990lbl 784 `"Solderers"', add
label define occ1990lbl 785 `"Assemblers of electrical equipment"', add
label define occ1990lbl 789 `"Hand painting, coating, and decorating occupations"', add
label define occ1990lbl 796 `"Production checkers and inspectors"', add
label define occ1990lbl 799 `"Graders and sorters in manufacturing"', add
label define occ1990lbl 803 `"Supervisors of motor vehicle transportation"', add
label define occ1990lbl 804 `"Truck, delivery, and tractor drivers"', add
label define occ1990lbl 808 `"Bus drivers"', add
label define occ1990lbl 809 `"Taxi cab drivers and chauffeurs"', add
label define occ1990lbl 813 `"Parking lot attendants"', add
label define occ1990lbl 823 `"Railroad conductors and yardmasters"', add
label define occ1990lbl 824 `"Locomotive operators (engineers and firemen)"', add
label define occ1990lbl 825 `"Railroad brake, coupler, and switch operators"', add
label define occ1990lbl 829 `"Ship crews and marine engineers"', add
label define occ1990lbl 834 `"Water transport infrastructure tenders and crossing guards"', add
label define occ1990lbl 844 `"Operating engineers of construction equipment"', add
label define occ1990lbl 848 `"Crane, derrick, winch, and hoist operators"', add
label define occ1990lbl 853 `"Excavating and loading machine operators"', add
label define occ1990lbl 859 `"Misc material moving occupations"', add
label define occ1990lbl 865 `"Helpers, constructions"', add
label define occ1990lbl 866 `"Helpers, surveyors"', add
label define occ1990lbl 869 `"Construction laborers"', add
label define occ1990lbl 874 `"Production helpers"', add
label define occ1990lbl 875 `"Garbage and recyclable material collectors"', add
label define occ1990lbl 876 `"Materials movers: stevedores and longshore workers"', add
label define occ1990lbl 877 `"Stock handlers"', add
label define occ1990lbl 878 `"Machine feeders and offbearers"', add
label define occ1990lbl 883 `"Freight, stock, and materials handlers"', add
label define occ1990lbl 885 `"Garage and service station related occupations"', add
label define occ1990lbl 887 `"Vehicle washers and equipment cleaners"', add
label define occ1990lbl 888 `"Packers and packagers by hand"', add
label define occ1990lbl 889 `"Laborers outside construction"', add
label define occ1990lbl 905 `"Military"', add
label define occ1990lbl 991 `"Unemployed"', add
label define occ1990lbl 999 `"Unknown"', add
label values occ1990 occ1990lbl

label define wkswork2lbl 0 `"N/A"'
label define wkswork2lbl 1 `"1-13 weeks"', add
label define wkswork2lbl 2 `"14-26 weeks"', add
label define wkswork2lbl 3 `"27-39 weeks"', add
label define wkswork2lbl 4 `"40-47 weeks"', add
label define wkswork2lbl 5 `"48-49 weeks"', add
label define wkswork2lbl 6 `"50-52 weeks"', add
label values wkswork2 wkswork2lbl

label define hrswork2lbl 0 `"N/A"'
label define hrswork2lbl 1 `"1-14 hours"', add
label define hrswork2lbl 2 `"15-29 hours"', add
label define hrswork2lbl 3 `"30-34 hours"', add
label define hrswork2lbl 4 `"35-39 hours"', add
label define hrswork2lbl 5 `"40 hours"', add
label define hrswork2lbl 6 `"41-48 hours"', add
label define hrswork2lbl 7 `"49-59 hours"', add
label define hrswork2lbl 8 `"60+ hours"', add
label values hrswork2 hrswork2lbl

save ipums_1990_main, replace


/***************************************************************************************************/
/***************************************************************************************************/
/************************  Step 2:  Create Main Data Files For Each Year  **************************/
/***************************************************************************************************/
/***************************************************************************************************/

*  Notes
*  For each year, we will make one analysis file that creates the key variables for our analysis.
*  These variables include employment status, occupation code, annual earnings, hourly wage
*  We also define the race*sex groups
*  The occupation code creation is the key part of this file.  We discuss these codes as we introduce them.
*  All files start in with raw Census or ACS files (as read in above)


/*  Step A:  Read in the File */    

	#delimit ;
	clear;
	clear matrix ;
	set mem 3000m;
	cd c:\ErikMain\discrimination_growth;
	use ipums_1990_main, clear; 
	

/*  Step B:  Restrict Sample */													

	keep if age >= 25 & age <= 54;   											/* Restrict to prime age individuals */
	keep if race == 1 | race == 2 ;  											/* Keep only black and white respondents */
	drop if occ1990 == 905 | empstatd == 13 | empstatd == 14 | empstatd == 15 ; /* Exclude the military */
	drop if occ1990 >= 999 & empstat == 1 ;   									/* Exclude those with missing occupation codes */ 
	drop if (empstatd == 20 | empstatd == 21 | empstatd == 22 | occ1990 == 991) /* Exclude the unemployed */; 


/* Step C:  Define Key Demographic Variables (employment status, groups, region, education, usual hours worked, etc.) */  

	
	/* C1.  Create work Status Codes */ ; 
													
	#delimit ;

	gen emp_full = empstat == 1 & uhrswork >= 30 ;					/* This is a measure of whether working full time */
	gen emp_part = empstat == 1 & uhrswork < 30 & uhrswork >= 15; 	/* This is a measure of whether working part time - we assume working less than 15 hours a week is not working */
	gen home = emp_full ~= 1 & emp_part ~= 1   ;					/* This is a measure of working in the home sector */

	gen 		emp_full_adj = . ;									/*emp_full_adj will be the variable we use for employed counts - it counts part time workers as 0.5 workers*/ 
	replace 	emp_full_adj = 0 	if home == 1 ;
	replace 	emp_full_adj = 0.5 	if emp_part == 1;
	replace 	emp_full_adj = 1 	if emp_full == 1;

	gen 		home_adj = . ; 										/* "home_adj" will be the variable use for home sector counts - it counts part time workers as being 0.5 in home sector */ 
	replace 	home_adj = 1 	if home == 1 ;
	replace 	home_adj = 0.5 	if emp_part == 1;
	replace 	home_adj = 0 	if emp_full == 1;

	gen 		person_adj = 1;
	replace 	person_adj = 0.5 if emp_part == 1 ; 

	gen emp_full_lastyear = wkswork1 >= 48 & (incwage + incbus + incfarm) >= 1000 * 129.9/207.16 ;   /*  This is to condition our earnings variables - it restricts sample to those working at least 48 weeks in prior year  */ 


	/*  C2.	Create Sex and Gender Controls */ ;

	gen white_man = race == 1 & sex == 1;
	gen white_woman = race == 1 & sex == 2;
	gen black_woman = race == 2 & sex == 2 ;
	gen black_man = race == 2 & sex == 1 ; 


	/*  C3.	Create Harmonized Education Controls */ ;

	gen highgrade = 0 if educ == 0 ;
	replace highgrade = 4 if educ == 1 ;
	replace highgrade = 8 if educ == 2 ;
	replace highgrade = 9 if educ == 3 ; 
	replace highgrade = 10 if educ == 4 ; 
	replace highgrade = 11 if educ == 5 ; 
	replace highgrade = 12 if educ == 6 ; 
	replace highgrade = 13 if educ == 7 ; 
	replace highgrade = 14 if educ == 8 ; 
	replace highgrade = 15 if educ == 9 ; 
	replace highgrade = 16 if educ == 10 ; 
	replace highgrade = 19 if educ == 11 ; 


/*  Step D:		Create Harmonized Wage/Earnings Data */ ; 

	gen incwage_full = .;     																			/*  This is our measure of annual earnings */ 	
	replace incwage_full = incwage + incbus + incfarm if emp_full_lastyear == 1 & emp_full_adj > 0;    	/*  Restrict this measure to only those who worked full time last year */           

	gen wage = incwage_full/(uhrswork * wkswork1) ;  	/*  We will use this measure to compute wage gaps between groups */ 

	gen ln_incwage_full = ln(incwage_full);
	gen ln_wage = ln(wage) ;

	gen incwage_full_adj = incwage_full ;        		/*  We will use this measure to compute occupational earnings */ 


/*  Step E - Create Occupation Codes  */ ; 

	*  Notes:   
	*  Occupation codes based on the 1990 Census Sub Categories of their occupation categories.  
	*  See http://usa.ipums.org/usa/volii/99occup.shtml for details  ;  

	/*   E1:  Create occ_code - this is our main harmonized set of occupation codes */ ;  
	
	gen occ_code = 0  ;  						
	replace occ_code = 1  if ((occ1990 >= 3 & occ1990 <= 22) ) & emp_full_adj > 0;
	replace occ_code = 2  if ((occ1990 >= 23 & occ1990 <= 37) | (occ1990 >= 303 & occ1990 <= 307) | (occ1990 == 200))& emp_full_adj > 0;
	replace occ_code = 3  if (occ1990 == 43)& emp_full_adj > 0;
	replace occ_code = 4  if ((occ1990 >= 44 & occ1990 <= 63) | occ1990 == 867)& emp_full_adj > 0;
	replace occ_code = 5  if (occ1990 >= 64 & occ1990 <= 68)& emp_full_adj > 0 ;
	replace occ_code = 6  if (occ1990 >= 69 & occ1990 <= 83)& emp_full_adj > 0;
	replace occ_code = 7  if (occ1990 >= 84 & occ1990 <= 89)& emp_full_adj > 0;
	replace occ_code = 8  if (occ1990 >= 95 & occ1990 <= 97)& emp_full_adj > 0;
	replace occ_code = 9  if (occ1990 >= 98 & occ1990 <= 106)& emp_full_adj > 0;
	replace occ_code = 10 if (occ1990 >= 113 & occ1990 <= 154)& emp_full_adj > 0;
	replace occ_code = 11  if (occ1990 >= 155 & occ1990 <= 163)& emp_full_adj > 0;
	replace occ_code = 12  if (occ1990 >= 164 & occ1990 <= 165)& emp_full_adj > 0;
	replace occ_code = 13  if (occ1990 >= 166 & occ1990 <= 173)& emp_full_adj > 0;
	replace occ_code = 14  if (occ1990 >= 174 & occ1990 <= 177)& emp_full_adj > 0;
	replace occ_code = 15  if (occ1990 >= 178 & occ1990 <= 179)& emp_full_adj > 0;
	replace occ_code = 16  if (occ1990 >= 183 & occ1990 <= 199)& emp_full_adj > 0;
	replace occ_code = 17  if (occ1990 >= 203 & occ1990 <= 208)& emp_full_adj > 0;
	replace occ_code = 18  if (occ1990 >= 213 & occ1990 <= 218)& emp_full_adj > 0;
	replace occ_code = 19  if (occ1990 >= 223 & occ1990 <= 225)& emp_full_adj > 0;
	replace occ_code = 20  if (occ1990 >= 226 & occ1990 <= 235)& emp_full_adj > 0;
	replace occ_code = 21  if (occ1990 >= 243 & occ1990 <= 290)& emp_full_adj > 0;
	replace occ_code = 22  if (occ1990 >= 313 & occ1990 <= 315)& emp_full_adj > 0;
	replace occ_code = 23  if (occ1990 >= 316 & occ1990 <= 323)& emp_full_adj > 0;
	replace occ_code = 24  if (occ1990 >= 325 & occ1990 <= 336)& emp_full_adj > 0;
	replace occ_code = 25  if (occ1990 >= 337 & occ1990 <= 344)& emp_full_adj > 0;
	replace occ_code = 26  if (occ1990 >= 345 & occ1990 <= 347)& emp_full_adj > 0;
	replace occ_code = 27  if ((occ1990 >= 348 & occ1990 <= 353) | (occ1990 >= 308 & occ1990 <= 309)) & emp_full_adj > 0;
	replace occ_code = 28  if (occ1990 >= 354 & occ1990 <= 357)& emp_full_adj > 0;
	replace occ_code = 29  if (occ1990 >= 359 & occ1990 <= 374)& emp_full_adj > 0;
	replace occ_code = 30  if (occ1990 >= 375 & occ1990 <= 378)& emp_full_adj > 0;
	replace occ_code = 31  if (occ1990 >= 379 & occ1990 <= 391)& emp_full_adj > 0;
	replace occ_code = 32  if (occ1990 >= 403 & occ1990 <= 408)& emp_full_adj > 0;
	replace occ_code = 33  if ((occ1990 >= 416 & occ1990 <= 417) | occ1990 == 413)& emp_full_adj > 0;   
	replace occ_code = 34  if ((occ1990 >= 418 & occ1990 <= 424) | occ1990 == 414)& emp_full_adj > 0;    
	replace occ_code = 35  if ((occ1990 >= 425 & occ1990 <= 427) | occ1990 == 415)& emp_full_adj > 0;     
	replace occ_code = 36  if (occ1990 >= 433 & occ1990 <= 444)& emp_full_adj > 0;
	replace occ_code = 37  if (occ1990 >= 445 & occ1990 <= 447)& emp_full_adj > 0;
	replace occ_code = 38  if (occ1990 >= 448 & occ1990 <= 455)& emp_full_adj > 0;
	replace occ_code = 39  if (occ1990 >= 456 & occ1990 <= 469)& emp_full_adj > 0;
	replace occ_code = 40  if (occ1990 >= 473 & occ1990 <= 476)& emp_full_adj > 0;
	replace occ_code = 41  if (occ1990 >= 477 & occ1990 <= 484)& emp_full_adj > 0;
	replace occ_code = 42  if (occ1990 >= 485 & occ1990 <= 489)& emp_full_adj > 0;
	replace occ_code = 43  if (occ1990 >= 494 & occ1990 <= 499)& emp_full_adj > 0;
	replace occ_code = 44  if (occ1990 >= 503 & occ1990 <= 519)& emp_full_adj > 0;
	replace occ_code = 45  if (occ1990 >= 523 & occ1990 <= 534)& emp_full_adj > 0;
	replace occ_code = 46  if (occ1990 >= 535 & occ1990 <= 549)& emp_full_adj > 0;
	replace occ_code = 47  if ((occ1990 >= 553 & occ1990 <= 599) | occ1990 == 866 | occ1990 == 869)& emp_full_adj > 0;
	replace occ_code = 48  if ((occ1990 >= 613 & occ1990 <= 617) | occ1990 == 868)& emp_full_adj > 0;
	replace occ_code = 49  if (occ1990 == 628)& emp_full_adj > 0;
	replace occ_code = 50  if (occ1990 >= 634 & occ1990 <= 655)& emp_full_adj > 0;
	replace occ_code = 51  if (occ1990 >= 656 & occ1990 <= 659)& emp_full_adj > 0;
	replace occ_code = 52  if (occ1990 >= 666 & occ1990 <= 674)& emp_full_adj > 0;
	replace occ_code = 53  if (occ1990 >= 675 & occ1990 <= 684)& emp_full_adj > 0;
	replace occ_code = 54  if (occ1990 >= 686 & occ1990 <= 688)& emp_full_adj > 0;
	replace occ_code = 55  if (occ1990 >= 694 & occ1990 <= 699)& emp_full_adj > 0;
	replace occ_code = 56  if (occ1990 >= 703 & occ1990 <= 717)& emp_full_adj > 0;
	replace occ_code = 57  if (occ1990 >= 719 & occ1990 <= 725)& emp_full_adj > 0;
	replace occ_code = 58  if (occ1990 >= 726 & occ1990 <= 733)& emp_full_adj > 0;
	replace occ_code = 59  if (occ1990 >= 734 & occ1990 <= 737)& emp_full_adj > 0;
	replace occ_code = 60  if (occ1990 >= 738 & occ1990 <= 749)& emp_full_adj > 0;
	replace occ_code = 61  if (occ1990 >= 753 & occ1990 <= 779)& emp_full_adj > 0;
	replace occ_code = 62  if ((occ1990 >= 783 & occ1990 <= 795) | occ1990 == 874) & emp_full_adj > 0 ;
	replace occ_code = 63  if (occ1990 >= 796 & occ1990 <= 799 | occ1990 >= 689 & occ1990 <= 693)& emp_full_adj > 0 ;
	replace occ_code = 64  if (occ1990 >= 803 & occ1990 <= 815)& emp_full_adj > 0;
	replace occ_code = 65  if ((occ1990 >= 823 & occ1990 <= 834) |(occ1990 >= 843 & occ1990 <= 865)) & emp_full_adj > 0;
	replace occ_code = 66  if (occ1990 >= 875 & occ1990 <= 890)& emp_full_adj > 0;


	label   define  occ_codelbl 0   `"Home"', add;
	label	define	occ_codelbl	1	`"Executives, Administrative, and Managerial"',	add;	
	label	define	occ_codelbl	2	`"Management Related"',	add;			
	label	define	occ_codelbl	3	`"Architects"',	add;				
	label	define	occ_codelbl	4	`"Engineers"',	add;				
	label	define	occ_codelbl	5	`"Math	and	Computer Science"',	add;	
	label	define	occ_codelbl	6	`"Natural Science"',	add;			
	label	define	occ_codelbl	7	`"Health Diagnosing"',	add;			
	label	define	occ_codelbl	8	`"Health Assessment"',	add;			
	label	define	occ_codelbl	9	`"Therapists"',	add;				
	label	define	occ_codelbl	10	`"Teachers,	Postsecondary"',	add;			
	label	define	occ_codelbl	11	`"Teachers,	Non-Postsecondary"',	add;			
	label	define	occ_codelbl	12	`"Librarians and Curators"',	add;		
	label	define	occ_codelbl	13	`"Social Scientists	and	Urban Planners"',	add;
	label	define	occ_codelbl	14	`"Social, Recreation, and Religious	Workers"',	add;
	label	define	occ_codelbl	15	`"Lawyers and	Judges"',	add;		
	label	define	occ_codelbl	16	`"Arts	and	Athletes"',	add;		
	label	define	occ_codelbl	17	`"Health Technicians"',	add;			
	label	define	occ_codelbl	18	`"Engineering Technicians"',	add;			
	label	define	occ_codelbl	19	`"Science Technicians"',	add;			
	label	define	occ_codelbl	20	`"Technicians, Other"',	add;			
	label	define	occ_codelbl	21	`"Sales, all"',	add;				
	label	define	occ_codelbl	22	`"Secretaries"',	add;				
	label	define	occ_codelbl	23	`"Information Clerks"',	add;			
	label	define	occ_codelbl	24	`"Records Processing, Non-Financial"',	add;		
	label	define	occ_codelbl	25	`"Records Processing, Financial"',	add;		
	label	define	occ_codelbl	26	`"Office Machine Operator"',	add;		
	label	define	occ_codelbl	27	`"Computer and Communications Equipment Operator"',	add;		
	label	define	occ_codelbl	28	`"Mail Distribution"',	add;			
	label	define	occ_codelbl	29	`"Scheduling and Distributing Clerks"',	add;	
	label	define	occ_codelbl	30	`"Adjusters	and	Investigators"',	add;		
	label	define	occ_codelbl	31	`"Misc.	Admin Support"',	add;		
	label	define	occ_codelbl	32	`"Private Household	Occupations"',	add;		
	label	define	occ_codelbl	33	`"Firefighting"',	add;				
	label	define	occ_codelbl	34	`"Police"',	add;				
	label	define	occ_codelbl	35	`"Guards"',	add;				
	label	define	occ_codelbl	36	`"Food Prep and Service"',	add;	
	label	define	occ_codelbl	37	`"Health Service"',	add;			
	label	define	occ_codelbl	38	`"Cleaning and Building	Service"',	add;	
	label	define	occ_codelbl	39	`"Personal Service"',	add;			
	label	define	occ_codelbl	40	`"Farm Managers"',	add;			
	label	define	occ_codelbl	41	`"Farm Non-Managers"',	add;			
	label	define	occ_codelbl	42	`"Related Agriculture"',	add;			
	label	define	occ_codelbl	43	`"Forest, Logging, Fishers and Hunter"',	add;				
	label	define	occ_codelbl	44	`"Vehicle Mechanic"',	add;			
	label	define	occ_codelbl	45	`"Electronic Repairer"',	add;			
	label	define	occ_codelbl	46	`"Misc.	Repairer"',	add;			
	label	define	occ_codelbl	47	`"Construction Trade"',	add;			
	label	define	occ_codelbl	48	`"Extractive"',	add;				
	label	define	occ_codelbl	49	`"Precision	Production,	Supervisor"',	add;		
	label	define	occ_codelbl	50	`"Precision	Metal"',	add;			
	label	define	occ_codelbl	51	`"Precision	Wood"',	add;			
	label	define	occ_codelbl	52	`"Precision, Textile"',	add;			
	label	define	occ_codelbl	53	`"Precision, Other"',	add;			
	label	define	occ_codelbl	54	`"Precision, Food"',	add;			
	label	define	occ_codelbl	55	`"Plant	and	System Operator"',	add;	
	label	define	occ_codelbl	56	`"Metal	and	Plastic	Machine	Operator"',	add;
	label	define	occ_codelbl	57	`"Metal	and	Plastic	Processing	Operator"',	add;
	label	define	occ_codelbl	58	`"Woodworking Machine Operator"',	add;		
	label	define	occ_codelbl	59	`"Textile Machine Operator"',	add;		
	label	define	occ_codelbl	60	`"Printing Machine Operator"',	add;		
	label	define	occ_codelbl	61	`"Machine Operator, Other"',	add;		
	label	define	occ_codelbl	62	`"Fabricators"',	add;				
	label	define	occ_codelbl	63	`"Production Inspectors"',	add;			
	label	define	occ_codelbl	64	`"Motor	Vehicle	Operator"',	add;		
	label	define	occ_codelbl	65	`"Non Motor Vehicle Operator"',	add;		
	label	define	occ_codelbl	66	`"Freight, Stock, Material	Handler"',	add;	

	label values occ_code occ_codelbl;
	
	
	/*	E2.  Create an alternate occupation code (used in creating data files for programs) */
	
	*	Notes:  The alternate occupation codes counts people working part time as being in the home sector (with 0.5 weight);
	
	gen occ_code2 = occ_code ;
	replace occ_code2 = 0 if person_adj == 0.5; 
	label values occ_code2 occ_codelbl ; 

	
	/*   E3: Create Broad Occupation Codes - These are only used for some results in the online appendix  */ ; 
	
	/* Note - Instead of defining 67 occupation - we only define 20 broad occupations */

	gen occ_broad = 0  ;  						
	replace occ_broad = 1  if ((occ1990 >= 3 & occ1990 <= 22) ) & emp_full_adj > 0;
	replace occ_broad = 2  if ((occ1990 >= 23 & occ1990 <= 37) | (occ1990 == 200))& emp_full_adj > 0;
	replace occ_broad = 3  if ((occ1990 >= 43 & occ1990 <= 68) | occ1990 == 867)& emp_full_adj > 0;
	replace occ_broad = 4  if ((occ1990 >= 69 & occ1990 <= 83) | (occ1990 >= 166 & occ1990 <= 177) | (occ1990 >= 183 & occ1990 <= 199)) & emp_full_adj > 0;
	replace occ_broad = 5  if ((occ1990 >= 84 & occ1990 <= 89) | (occ1990 >= 178 & occ1990 <= 179)) & emp_full_adj > 0 ;
	replace occ_broad = 6  if ((occ1990 >= 95 & occ1990 <= 106) | (occ1990 >= 445 & occ1990 <= 447))  & emp_full_adj > 0;
	replace occ_broad = 7  if (occ1990 >= 113 & occ1990 <= 154) & emp_full_adj > 0;
	replace occ_broad = 8  if (occ1990 >= 155 & occ1990 <= 165) & emp_full_adj > 0;
	replace occ_broad = 9  if (occ1990 >= 203 & occ1990 <= 235) & emp_full_adj > 0;
	replace occ_broad = 10 if (occ1990 >= 243 & occ1990 <= 290) & emp_full_adj > 0;
	replace occ_broad = 11 if (occ1990 >= 303 & occ1990 <= 391) & emp_full_adj > 0;
	replace occ_broad = 12 if (occ1990 >= 413 & occ1990 <= 427) & emp_full_adj > 0;
	replace occ_broad = 13 if ((occ1990 >= 403 & occ1990 <= 408) | (occ1990 >= 433 & occ1990 <= 444)| (occ1990 >= 448 & occ1990 <= 469)) & emp_full_adj > 0;
	replace occ_broad = 14  if ((occ1990 >= 473 & occ1990 <= 499) | (occ1990 >= 613 & occ1990 <= 617) | (occ1990 == 868)) & emp_full_adj > 0;
	replace occ_broad = 15  if ((occ1990 >= 503 & occ1990 <= 599)| (occ1990 == 866) | (occ1990 == 869)) & emp_full_adj > 0;
	replace occ_broad = 16  if (occ1990 >= 628 & occ1990 <= 688) & emp_full_adj > 0;
	replace occ_broad = 17  if (occ1990 >= 694 & occ1990 <= 779) & emp_full_adj > 0;
	replace occ_broad = 18  if ((occ1990 >= 783 & occ1990 <= 799) | (occ1990 >= 689 & occ1990 <= 693) |(occ1990 >= 875 & occ1990 <= 890) | (occ1990 == 874)) & emp_full_adj > 0;
	replace occ_broad = 19  if ((occ1990 >= 803 & occ1990 <= 815) | (occ1990 >= 823 & occ1990 <= 834) |(occ1990 >= 843 & occ1990 <= 865))& emp_full_adj > 0;

	label   define  occ_broadlbl 	0   `"Home"', add;
	label	define	occ_broadlbl	1	`"Executives, Administrative, and Managerial"',	add;	
	label	define	occ_broadlbl	2	`"Management Related"',	add;			
	label	define	occ_broadlbl	3	`"Architects, Engineers, Math, and Computer Science"',	add;				
	label	define	occ_broadlbl	4	`"Natural and Social Scientists, Recreation, Religious, Arts, Athletes"',	add;	
	label	define	occ_broadlbl	5	`"Doctors and Lawyers"',	add;				
	label	define	occ_broadlbl	6	`"Nurses, Therapists, and Other Health Service"',	add;	
	label	define	occ_broadlbl	7	`"Teachers, Postsecondary"',	add;			
	label	define	occ_broadlbl	8	`"Teachers, Non-Postsecondary and Librarians"',	add;			
	label	define	occ_broadlbl 	9	`"Health and Science Technicians"',	add;			
	label	define	occ_broadlbl	10	`"Sales, All"',	add;				
	label	define	occ_broadlbl	11	`"Administrative Support, Clerks, Record Keepers"',	add;			
	label	define	occ_broadlbl	12	`"Fire, Police, and Guards"',	add;			
	label	define	occ_broadlbl	13	`"Food, Cleaning, and Personal Services and Private Household"',	add;		
	label	define	occ_broadlbl	14	`"Farm, Related Agrigulture, Logging, and Extraction"',	add;
	label	define	occ_broadlbl	15	`"Mechanics and Construction"',	add;
	label	define	occ_broadlbl	16	`"Precision Manufacturing"',	add;		
	label	define	occ_broadlbl	17	`"Manufacturing Operators"',	add;		
	label	define	occ_broadlbl	18	`"Fabricators, Inspectors, and Material Handlers"',	add;			
	label	define	occ_broadlbl	19 `"Vehicle Operators"',	add;			
			

	label values occ_broad occ_broadlbl;

	gen occ_broad2 = occ_broad ;
	replace occ_broad2 = 0 if person_adj == 0.5; 

	label values occ_broad2 occ_broadlbl ; 

	
/*   Part F:  Save the extract as a new file */ ; 

	#delimit ; 
	save 1990_extract_composite_main, replace ; 





