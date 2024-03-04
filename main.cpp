/*
 * This file is part of the VanitySearch distribution (https://github.com/JeanLucPons/VanitySearch).
 * Copyright (c) 2019 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "Timer.h"
#include "Vanity.h"
#include "SECP256k1.h"
#include <fstream>
#include <string>
#include <string.h>
#include <stdexcept>
#include "hash/sha512.h"
#include "hash/sha256.h"

#define RELEASE "1.16 Linux with BitCrack integration"

using namespace std;

// ------------------------------------------------------------------------------------------

void printUsage() {

  printf("VanitySeacrh [-check] [-v] [-u] [-b] [-c] [-gpu] [-stop] [-i inputfile]\n");
  printf("             [-gpuId gpuId1[,gpuId2,...]] [-g gridSize1[,gridSize2,...]]\n");
  printf("             [-o outputfile] [-m maxFound] [-ps seed] [-s seed] [-t nbThread]\n");
  printf("             [-nosse] [-r rekey] [-check] [-kp] [-sp startPubKey]\n");
  printf("             [-rp privkey partialkeyfile] [prefix]\n\n");
  printf(" prefix: prefix to search (Can contains wildcard '?' or '*')\n");
  printf(" -v: Print version\n");
  printf(" -u: Search uncompressed addresses\n");
  printf(" -b: Search both uncompressed or compressed addresses\n");
  printf(" -c: Case unsensitive search\n");
  printf(" -gpu: Enable gpu calculation\n");
  printf(" -stop: Stop when all prefixes are found\n");
  printf(" -i inputfile: Get list of prefixes to search from specified file\n");
  printf(" -o outputfile: Output results to the specified file\n");
  printf(" -gpu gpuId1,gpuId2,...: List of GPU(s) to use, default is 0\n");
  printf(" -g gridSize1,gridSize2,...: Specify GPU(s) kernel gridsize, default is 8*(MP number)\n");
  printf(" -m: Specify maximun number of prefixes found by each kernel call\n");
  printf(" -t threadNumber: Specify number of CPU thread, default is number of core\n");
  printf(" -nosse: Disable SSE hash function\n");
  printf(" -l: List cuda enabled devices\n");
  printf(" -check: Check CPU and GPU kernel vs CPU\n");
  printf(" -kp: Generate key pair\n");
  printf(" -rp privkey partialkeyfile: Reconstruct final private key(s) from partial key(s) info.\n");
  printf(" -sp startPubKey: Start the search with a pubKey (for private key splitting)\n");
  printf(" -r rekey: Rekey interval in MegaKey, cpu=1000 gpu=100000\n");
  printf("\n");
  printf(" --keyspace START \n");
  printf("            START:END \n");
  printf("            START:+COUNT \n");
  printf("            :+COUNT \n");
  printf("            :END \n");
  printf("            Where START, END, COUNT are in hex format\n");
  printf(" --share M / N: Divide the keyspace into N equal shares, process the Mth share\n");
  printf(" --continue sessfile: Save / load progress from specified file\n");
  exit(-1);

}

// ------------------------------------------------------------------------------------------

int getInt(string name,char *v) {

  int r;

  try {

    r = std::stoi(string(v));

  } catch(std::invalid_argument&) {

    printf("[ERROR] Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

  return r;

}

// ------------------------------------------------------------------------------------------

void getInts(string name,vector<int> &tokens, const string &text, char sep) {

  size_t start = 0, end = 0;
  tokens.clear();
  int item;

  try {

    while ((end = text.find(sep, start)) != string::npos) {
      item = std::stoi(text.substr(start, end - start));
      tokens.push_back(item);
      start = end + 1;
    }

    item = std::stoi(text.substr(start));
    tokens.push_back(item);

  } catch(std::invalid_argument &) {

    printf("[ERROR] Invalid %s argument, number expected\n",name.c_str());
    exit(-1);

  }

}

// ------------------------------------------------------------------------------------------

void getKeySpace( const string &text, BITCRACK_PARAM * bc, Int& maxKey) {

	size_t start = 0, end = 0;
	string item;

	try {

		if((end = text.find(':', start)) != string::npos) {
			item = std::string(text.substr(start, end));
			start = end + 1;
		}
		else {
			item = std::string(text);
		}

		if (item.length() == 0) {
			bc->ksStart.SetInt32(1);
		}
		else if (item.length() > 64) {
			printf("[ERROR] keyspaceSTART: invalid privkey (64 length)\n");
			exit(-1);
		}
		else{
			item.insert(0, 64-item.length(), '0');
			for (int i = 0; i < 32; i++) {
				unsigned char my1ch = 0;
				sscanf(&item[2 * i], "%02hhX", &my1ch);
				bc->ksStart.SetByte(31 - i, my1ch);
			}
		}

		if (start != 0 && (end = text.find('+', start)) != string::npos) {
			item = std::string(text.substr(end + 1));
			if (item.length() > 64 || item.length() == 0) {
				printf("[ERROR] keyspace__END: invalid privkey (64 length)\n");
				exit(-1);
			}
			item.insert(0, 64 - item.length(), '0');
			for (int i = 0; i < 32; i++) {
				unsigned char my1ch = 0;
				sscanf(&item[2 * i], "%02hhX", &my1ch);
				bc->ksFinish.SetByte(31 - i, my1ch);
			}
			bc->ksFinish.Add(&bc->ksStart);
		}
		else if (start != 0) {
			item = std::string(text.substr(start));
			if (item.length() > 64 || item.length() == 0) {
				printf("[ERROR] keyspace__END: invalid privkey (64 length)\n");
				exit(-1);
			}
			item.insert(0, 64 - item.length(), '0');
			for (int i = 0; i < 32; i++) {
				unsigned char my1ch = 0;
				scanf(&item[2 * i], "%02hhX", &my1ch);
				bc->ksFinish.SetByte(31 - i, my1ch);
			}
		}
		else {
			bc->ksFinish.Set(&maxKey);
		}

	}
	catch (std::invalid_argument &) {

		printf("[ERROR] Invalid --keyspace argument \n");
		exit(-1);

	}

}
// ------------------------------------------------------------------------------------------

void checkKeySpace(BITCRACK_PARAM * bc, Int& maxKey) {

	if (bc->ksStart.IsGreater(&maxKey) || bc->ksFinish.IsGreater(&maxKey)) {
		printf("[ERROR] START/END IsGreater %064s \n", maxKey.GetBase16().c_str());
		exit(-1);
	}
	if (bc->ksFinish.IsLowerOrEqual(&bc->ksStart)) {
		printf("[ERROR] END IsLowerOrEqual START \n");
		exit(-1);
	}
	if (bc->ksFinish.IsLowerOrEqual(&bc->ksNext)) {
		printf("[ERROR] END: IsLowerOrEqual NEXT \n");
		exit(-1);
	}

	return;
}
// ------------------------------------------------------------------------------------------

void getShare(const string &text, int *shareM, int *shareN) {

	size_t start = 0, end = 0;
	string item;

	try {

		if ((end = text.find('/', start)) != string::npos) {
			*shareM = std::stoi(text.substr(start, end));
			start = end + 1;
			*shareN = std::stoi(text.substr(start));
		}
		else {
			printf("[ERROR] Invalid --share M/N argument \n");
			exit(-1);
		}
		if(*shareM <= 0 || *shareN <= 0 || *shareM > *shareN ) {
			printf("[ERROR] Invalid --share argument, need M<=N and M>0 and N>0 \n");
			exit(-1);
		}

	}
	catch (std::invalid_argument &) {

		printf("[ERROR] Invalid --share argument \n");
		exit(-1);

	}

}
// ------------------------------------------------------------------------------------------

void loadProgress(string fileName, BITCRACK_PARAM * bc) {

	if (fileName.length() > 0) {

		FILE *fp;
		fp = fopen(fileName.c_str(), "r");
		if (fp != NULL) {
			printf("[load] from sessfile '%s' \n", fileName.c_str());

			char f_buf[100];
			string f_str;
			size_t f_start = 0, f_end = 0;

			while (!feof(fp)) {
				if (fgets(f_buf, 99, fp)) {
					f_str = f_buf;
					if ((f_end = f_str.find("start=", f_start)) != string::npos) {
						f_str = std::string(f_str.substr(f_end + 6, 64));
						bc->ksStart.SetBase16((char *)f_str.c_str());
					}
					if ((f_end = f_str.find("next=", f_start)) != string::npos) {
						f_str = std::string(f_str.substr(f_end + 5, 64));
						bc->ksNext.SetBase16((char *)f_str.c_str());
					}
					if ((f_end = f_str.find("end=", f_start)) != string::npos) {
						f_str = std::string(f_str.substr(f_end + 4, 64));
						bc->ksFinish.SetBase16((char *)f_str.c_str());
					}
				}
			}
			printf("[load] start=%064s \n", bc->ksStart.GetBase16().c_str());
			printf("[load]  next=%064s \n", bc->ksNext.GetBase16().c_str());
			printf("[load]   end=%064s \n", bc->ksFinish.GetBase16().c_str());
			fclose(fp);
		}
	}

}
// ------------------------------------------------------------------------------------------

void parseFile(string fileName, vector<string> &lines) {

  // Get file size
  FILE *fp = fopen(fileName.c_str(), "rb");
  if (fp == NULL) {
    printf("[ERROR] ParseFile: cannot open %s %s\n", fileName.c_str(), strerror(errno));
    exit(-1);
  }
  fseek(fp, 0L, SEEK_END);
  size_t sz = ftell(fp);
  size_t nbAddr = sz / 33; /* Upper approximation */
  bool loaddingProgress = sz > 100000;
  fclose(fp);

  // Parse file
  int nbLine = 0;
  string line;
  ifstream inFile(fileName);
  lines.reserve(nbAddr);
  while (getline(inFile, line)) {

    // Remove ending \r\n
    int l = (int)line.length() - 1;
    while (l >= 0 && isspace(line.at(l))) {
      line.pop_back();
      l--;
    }

    if (line.length() > 0) {
      lines.push_back(line);
      nbLine++;
      if (loaddingProgress) {
        if ((nbLine % 50000) == 0)
          printf("[Loading input file %5.1f%%]\r", ((double)nbLine*100.0) / ((double)(nbAddr)*33.0 / 34.0));
      }
    }

  }

  if (loaddingProgress)
    printf("[Loading input file 100.0%%]\n");

}
// ------------------------------------------------------------------------------------------

void generateKeyPair(Secp256K1 *secp, string seed, int searchMode,bool paranoiacSeed) {

  if (seed.length() < 8) {
    printf("[ERROR] Use a seed of at leats 8 characters to generate a key pair\n");
    printf("Ex: VanitySearch -s \"A Strong Password\" -kp\n");
    exit(-1);
  }

  if(paranoiacSeed)
    seed = seed + Timer::getSeed(32);

  if (searchMode == SEARCH_BOTH) {
    printf("[ERROR] Use compressed or uncompressed to generate a key pair\n");
    exit(-1);
  }

  bool compressed = (searchMode == SEARCH_COMPRESSED);

  string salt = "0";
  unsigned char hseed[64];
  pbkdf2_hmac_sha512(hseed, 64, (const uint8_t *)seed.c_str(), seed.length(),
    (const uint8_t *)salt.c_str(), salt.length(),
    2048);

  Int privKey;
  privKey.SetInt32(0);
  sha256(hseed, 64, (unsigned char *)privKey.bits64);
  Point p = secp->ComputePublicKey(&privKey);
  printf("Priv : %s\n", secp->GetPrivAddress(compressed,privKey).c_str());
  printf("Pub  : %s\n", secp->GetPublicKeyHex(compressed,p).c_str());

}

// ------------------------------------------------------------------------------------------

void outputAdd(string outputFile, int addrType, string addr, string pAddr, string pAddrHex) {

  FILE *f = stdout;
  bool needToClose = false;

  if (outputFile.length() > 0) {
    f = fopen(outputFile.c_str(), "a");
    if (f == NULL) {
      printf("Cannot open %s for writing\n", outputFile.c_str());
      f = stdout;
    } else {
      needToClose = true;
    }
  }

  fprintf(f, "\nPub Addr: %s\n", addr.c_str());


  switch (addrType) {
  case P2PKH:
    fprintf(f, "Priv (WIF): p2pkh:%s\n", pAddr.c_str());
    break;
  case P2SH:
    fprintf(f, "Priv (WIF): p2wpkh-p2sh:%s\n", pAddr.c_str());
    break;
  case BECH32:
    fprintf(f, "Priv (WIF): p2wpkh:%s\n", pAddr.c_str());
    break;
  }
  fprintf(f, "Priv (HEX): 0x%064s\n", pAddrHex.c_str());

  if (needToClose)
    fclose(f);

}

// ------------------------------------------------------------------------------------------
#define CHECK_ADDR()                                           \
  fullPriv.ModAddK1order(&e, &partialPrivKey);                 \
  p = secp->ComputePublicKey(&fullPriv);                       \
  cAddr = secp->GetAddress(addrType, compressed, p);           \
  if (cAddr == addr) {                                         \
    found = true;                                              \
    string pAddr = secp->GetPrivAddress(compressed, fullPriv); \
    string pAddrHex = fullPriv.GetBase16();                    \
    outputAdd(outputFile, addrType, addr, pAddr, pAddrHex);    \
  }

void reconstructAdd(Secp256K1 *secp, string fileName, string outputFile, string privAddr) {

  bool compressed;
  int addrType;
  Int lambda;
  Int lambda2;
  lambda.SetBase16("5363ad4cc05c30e0a5261c028812645a122e22ea20816678df02967c1b23bd72");
  lambda2.SetBase16("ac9c52b33fa3cf1f5ad9e3fd77ed9ba4a880b9fc8ec739c2e0cfc810b51283ce");

  Int privKey = secp->DecodePrivateKey((char *)privAddr.c_str(),&compressed);
  if(privKey.IsNegative())
    exit(-1);

  vector<string> lines;
  parseFile(fileName,lines);

  for (int i = 0; i < (int)lines.size(); i+=2) {

    string addr;
    string partialPrivAddr;

    if (lines[i].substr(0, 10) == "Pub Addr: ") {

      addr = lines[i].substr(10);

      switch (addr.data()[0]) {
      case '1':
        addrType = P2PKH; break;
      case '3':
        addrType = P2SH; break;
      case 'b':
      case 'B':
        addrType = BECH32; break;
      default:
        printf("Invalid partialkey info file at line %d\n", i);
        printf("%s Address format not supported\n", addr.c_str());
        continue;
      }

    } else {
      printf("[ERROR] Invalid partialkey info file at line %d (\"Pub Addr: \" expected)\n",i);
      exit(-1);
    }

    if (lines[i+1].substr(0, 13) == "PartialPriv: ") {
      partialPrivAddr = lines[i+1].substr(13);
    } else {
      printf("[ERROR] Invalid partialkey info file at line %d (\"PartialPriv: \" expected)\n", i);
      exit(-1);
    }

    bool partialMode;
    Int partialPrivKey = secp->DecodePrivateKey((char *)partialPrivAddr.c_str(), &partialMode);
    if (privKey.IsNegative()) {
      printf("[ERROR] Invalid partialkey info file at line %d\n", i);
      exit(-1);
    }

    if (partialMode != compressed) {

      printf("[WARNING] Invalid partialkey at line %d (Wrong compression mode, ignoring key)\n", i);
      continue;

    } else {

      // Reconstruct the address
      Int fullPriv;
      Point p;
      Int e;
      string cAddr;
      bool found = false;

      // No sym, no endo
      e.Set(&privKey);
      CHECK_ADDR();


      // No sym, endo 1
      e.Set(&privKey);
      e.ModMulK1order(&lambda);
      CHECK_ADDR();

      // No sym, endo 2
      e.Set(&privKey);
      e.ModMulK1order(&lambda2);
      CHECK_ADDR();

      // sym, no endo
      e.Set(&privKey);
      e.Neg();
      e.Add(&secp->order);
      CHECK_ADDR();

      // sym, endo 1
      e.Set(&privKey);
      e.ModMulK1order(&lambda);
      e.Neg();
      e.Add(&secp->order);
      CHECK_ADDR();

      // sym, endo 2
      e.Set(&privKey);
      e.ModMulK1order(&lambda2);
      e.Neg();
      e.Add(&secp->order);
      CHECK_ADDR();

      if (!found) {
        printf("Unable to reconstruct final key from partialkey line %d\n Addr: %s\n PartKey: %s\n",
          i, addr.c_str(),partialPrivAddr.c_str());
      }

    }

  }

}

// ------------------------------------------------------------------------------------------

int main(int argc, char* argv[]) {

  // Global Init
  Timer::Init();
  rseed((unsigned long)time(NULL));

  // Init SecpK1
  Secp256K1 *secp = new Secp256K1();
  secp->Init();

  // Browse arguments
  if (argc < 2) {
    printf("Not enough argument\n");
    printUsage();
  }

  int a = 1;
  bool gpuEnable = false;
  bool stop = false;
  int searchMode = SEARCH_COMPRESSED;
  vector<int> gpuId = {0};
  vector<int> gridSize = {-1};
  string seed = "";
  vector<string> prefix;
  string outputFile = "";
  int nbCPUThread = Timer::getCoreNumber();
  bool tSpecified = false;
  bool sse = true;
  uint32_t maxFound = 65536;
  uint64_t rekey = 1000;
  Point startPuKey;
  startPuKey.Clear();
  bool startPubKeyCompressed;
  bool caseSensitive = true;
  bool paranoiacSeed = false;

  string sessFile = "";

  Int maxKey;
  maxKey.SetBase16("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364140");

  BITCRACK_PARAM bitcrack, *bc;
  bc = &bitcrack;
  bc->ksStart.SetInt32(1);
  bc->ksNext.Set(&bc->ksStart);
  bc->ksFinish.Set(&maxKey);
  bc->shareM = 1;
  bc->shareN = 1;


  while (a < argc) {

    if (strcmp(argv[a], "-gpu")==0) {
      gpuEnable = true;
      a++;
    } else if (strcmp(argv[a], "-gpuId")==0) {
      a++;
      getInts("gpuId",gpuId,string(argv[a]),',');
      a++;
    } else if (strcmp(argv[a], "-stop") == 0) {
      stop = true;
      a++;
    } else if (strcmp(argv[a], "-c") == 0) {
      caseSensitive = false;
      a++;
    } else if (strcmp(argv[a], "-v") == 0) {
      printf("%s\n",RELEASE);
      exit(0);
    } else if (strcmp(argv[a], "-check") == 0) {

      Int::Check();
      secp->Check();

#ifdef WITHGPU
      GPUEngine g(gridSize[0],gpuId[0],maxFound,false);
      g.SetSearchMode(searchMode);
      g.Check(secp);
#else
  printf("GPU code not compiled, use -DWITHGPU when compiling.\n");
#endif
      exit(0);
    } else if (strcmp(argv[a], "-l") == 0) {

#ifdef WITHGPU
      GPUEngine::PrintCudaInfo();
#else
  printf("GPU code not compiled, use -DWITHGPU when compiling.\n");
#endif
      exit(0);

    } else if (strcmp(argv[a], "-kp") == 0) {
      generateKeyPair(secp,seed,searchMode,paranoiacSeed);
      exit(0);
    } else if (strcmp(argv[a], "-sp") == 0) {
      a++;
      string pub = string(argv[a]);
      startPuKey = secp->ParsePublicKeyHex(pub, startPubKeyCompressed);
      a++;
    } else if (strcmp(argv[a], "-rp") == 0) {
      a++;
      string priv = string(argv[a]);
      a++;
      string file = string(argv[a]);
      a++;
      reconstructAdd(secp,file,outputFile,priv);
      exit(0);
    } else if (strcmp(argv[a], "-u") == 0) {
      searchMode = SEARCH_UNCOMPRESSED;
      a++;
    } else if (strcmp(argv[a], "-b") == 0) {
      searchMode = SEARCH_BOTH;
      a++;
    } else if (strcmp(argv[a], "-nosse") == 0) {
      sse = false;
      a++;
    } else if (strcmp(argv[a], "-g") == 0) {
      a++;
      getInts("gridSize",gridSize,string(argv[a]),',');
      a++;
    } else if (strcmp(argv[a], "-s") == 0) {
      a++;
      seed = string(argv[a]);
      a++;
    } else if (strcmp(argv[a], "-ps") == 0) {
      a++;
      seed = string(argv[a]);
      paranoiacSeed = true;
      a++;
    } else if (strcmp(argv[a], "-o") == 0) {
      a++;
      outputFile = string(argv[a]);
      a++;
    } else if (strcmp(argv[a], "-i") == 0) {
      a++;
      parseFile(string(argv[a]),prefix);
      a++;
    } else if (strcmp(argv[a], "-t") == 0) {
      a++;
      nbCPUThread = getInt("nbCPUThread",argv[a]);
      a++;
      tSpecified = true;
    } else if (strcmp(argv[a], "-m") == 0) {
      a++;
      maxFound = getInt("maxFound", argv[a]);
      a++;
    } else if (strcmp(argv[a], "-r") == 0) {
      a++;
      rekey = (uint64_t)getInt("rekey", argv[a]);
      a++;
	}
	else if (strcmp(argv[a], "--continue") == 0) {
		a++;
		sessFile = string(argv[a]);
		a++;
	}
	else if (strcmp(argv[a], "--keyspace") == 0) {
		a++;
		getKeySpace(string(argv[a]), bc, maxKey);
		bc->ksNext.Set(&bc->ksStart);
		checkKeySpace(bc, maxKey);
		a++;
	}
	else if (strcmp(argv[a], "--share") == 0) {
		a++;
		getShare(string(argv[a]), &bc->shareM, &bc->shareN);
		a++;
	}
	else if (strcmp(argv[a], "--generator") == 0) {
		a++;
		int newgrpsize = 1024;
		newgrpsize = getInt("newgrpsize", argv[a]);
		GPUEngine::GenerateCode(secp, newgrpsize);
		printf("\n[Created GPUGroup.h][NEW_GRP_SIZE=%i]\n",newgrpsize);
		exit(0);
		a++;
    } else if (a == argc - 1) {
      prefix.push_back(string(argv[a]));
      a++;
    } else {
      printf("Unexpected %s argument\n",argv[a]);
      printUsage();
    }

  }

  printf("VanitySearch v" RELEASE "\n");

  if(gpuId.size()!=gridSize.size()) {
    if(gridSize.size()==1 && gridSize[0]==-1) {
      gridSize.clear();
      for(int i=0;i<gpuId.size();i++)
        gridSize.push_back(-1);
    } else {
      printf("Invalid gridSize or gpuId argument, must have same size\n");
      printUsage();
    }
  }

  // Let one CPU core free per gpu is gpu is enabled
  // It will avoid to hang the system
  if( !tSpecified && nbCPUThread>1 && gpuEnable)
    nbCPUThread-=(int)gpuId.size();
  if(nbCPUThread<0)
    nbCPUThread = 0;

  if (rekey == 0) rekey = 1000;
  if (gpuEnable && rekey <= 100000) rekey = 100000;


  // If a starting public key is specified, force the search mode according to the key
  if (!startPuKey.isZero()) {
	  searchMode = (startPubKeyCompressed) ? SEARCH_COMPRESSED : SEARCH_UNCOMPRESSED;
  }


  //Share to keyspace - apply before load progress

  if(bc->shareN > 1){
	  //printf("[share][before] start=%064s\n", bc->ksStart.GetBase16().c_str());
	  //printf("[share][before]   end=%064s\n", bc->ksFinish.GetBase16().c_str());

	  Int shareRange;
	  shareRange.Sub(&bc->ksFinish,&bc->ksStart);
	  Int BNshareN;
	  BNshareN.SetInt32((uint32_t)bc->shareN);
	  Int share1Key;
	  share1Key.Set(&shareRange);
	  share1Key.Div(&BNshareN);
	  share1Key.AddOne();
	  for (int i = 1; i < bc->shareM; i++) {
		  bc->ksStart.Add(&share1Key);
	  }
	  bc->ksFinish.Set(&bc->ksStart);
	  bc->ksFinish.Add(&share1Key);

	  bc->ksNext.Set(&bc->ksStart);

	  printf("[share] %i/%i \n", bc->shareM, bc->shareN);

  }


  // Try load progress from sessfile
  loadProgress(sessFile, bc);
  checkKeySpace(bc, maxKey);

  printf("[keyspace] start=%064s\n", bc->ksStart.GetBase16().c_str());
  printf("[keyspace]   end=%064s\n", bc->ksFinish.GetBase16().c_str());


  VanitySearch *v = new VanitySearch(secp, prefix, seed, searchMode, gpuEnable, stop, outputFile, sse,
    maxFound, rekey, caseSensitive, startPuKey, paranoiacSeed, sessFile, bc);
  v->Search(nbCPUThread,gpuId,gridSize);


  return 0;
}
