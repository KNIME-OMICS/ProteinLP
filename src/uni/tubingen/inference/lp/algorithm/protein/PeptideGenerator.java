package uni.tubingen.inference.lp.algorithm.protein;

import java.io.*;
import java.util.*;

/**
 * Date: Sep 14, 2007
 * Time: 22:11:14 PM
 *
 */
public class PeptideGenerator implements Runnable
{
    public static double[] AMINO_ACID_AVERAGE_MASSES = getMasses(false);
    public static double[] AMINO_ACID_MONOISOTOPIC_MASSES = getMasses(true);
    //Index of H-ION in acid tables
    public static final int H_ION_INDEX = 0;
    public static final double ELECTRON_MASS = 5.485e-4;

    public static final int DIGEST_TRYPTIC = 1;
    public static final int DIGEST_ALL = 2;

    private String m_inputFileName;
    private String m_outputFileName;
    private boolean m_computePI = true;
    private boolean m_computeHp = true;
    private int m_hpWindowSize = 9;
    private boolean m_computeAverageMass = true;
    private boolean m_computeMonoisotopicMass = true;
    private int m_digest=DIGEST_TRYPTIC;
    private boolean m_countOnly;
    private int m_maxMissedCleavages = 0;
    private double m_minMass = 1;
    private double m_maxMass = 640000;
    private int m_maxProteins = Integer.MAX_VALUE;
    private int m_minResidues = -1;
    private int m_maxResidues = Integer.MAX_VALUE;
    private double[] m_massTab = AMINO_ACID_MONOISOTOPIC_MASSES;
    private boolean m_async = true;
    private ArrayList m_peptides;

    private ArrayList m_listeners = new ArrayList();

    private int m_protNum = 0;
    private long m_pepCount = 0;
    private long m_aaCount = 0;

    /*You can test the functions via this main function*/
    public static void main(String[] args)
    {
        PeptideGenerator gp = initFromParams(args);
        if (null != gp)
        {
            if (!gp.m_countOnly)
            {
                System.out.print("prot");
                System.out.print('\t');

                System.out.print("pep");

                if (gp.m_computeAverageMass)
                {
                    System.out.print('\t');
                    System.out.print("m");
                }
                if (gp.m_computeMonoisotopicMass)
                {
                    System.out.print('\t');
                    System.out.print("mm");
                }
                if (gp.m_computePI)
                {
                    System.out.print('\t');
                    System.out.print("pi");
                }
                if(gp.m_computeHp)
                {
                    System.out.print('\t');
                    System.out.print("hp");
                }
                System.out.println();
            }

            Date d = new Date();
            gp.run();

            if (gp.m_countOnly)
            {
                System.out.print("Proteins: ");
                System.out.println(gp.m_protNum);
                System.out.print("Peptides: ");
                System.out.println(gp.m_pepCount);
                System.out.print("Residues: ");
                System.out.println(gp.m_aaCount);
            }
            System.err.println("Time: " + String.valueOf(new Date().getTime() - d.getTime()));

        }
    }

    public static PeptideGenerator initFromParams(String[] params)
    {
        if (params.length == 0)
            return paramError();

        PeptideGenerator gp = new PeptideGenerator();
        for (int i = 0; i < params.length; i++)
        {
            if (params[i].startsWith("-compute"))
            {
                HashSet computeSet = new HashSet();
                String[] split1 = params[i].split("=");
                if (split1.length != 2)
                    return paramError();

                String[] split2 = split1[1].split(",");
                for (int j = 0; j < split2.length; j++)
                    computeSet.add(split2[j].trim());

                gp.m_computeAverageMass =  computeSet.contains("m");
                gp.m_computeMonoisotopicMass = computeSet.contains("mm");
                gp.m_computeHp = computeSet.contains("hp");
                gp.m_computePI = computeSet.contains("pi");
            }
            else if (params[i].startsWith("-digest"))
            {
                String[] split = params[i].split("=");
                if (split.length != 2)
                    return paramError();
                String digest = split[1];
                if ("tryptic".equals(digest))
                    gp.m_digest = DIGEST_TRYPTIC;
                else if ("all".equals(digest))
                    gp.m_digest = DIGEST_ALL;
                else
                    return paramError();
            }
            else if (params[i].startsWith("-countOnly"))
            {
                gp.m_countOnly = true;
            }
            else if (params[i].startsWith("-maxProteins"))
            {
                String[] split = params[i].split("=");
                if (split.length != 2)
                    return paramError();
                gp.m_maxProteins = Integer.parseInt(split[1]);
            }
            else if (params[i].startsWith("-massRange"))
            {
                String[] split = params[i].split("=");
                if (split.length != 2)
                    return paramError();
                String[] massStrings = split[1].split("-");
                gp.m_minMass = Double.parseDouble(massStrings[0]);
                gp.m_maxMass = Double.parseDouble(massStrings[1]);
            }
            else if (params[i].startsWith("-maxMissed"))
            {
                String[] split = params[i].split("=");
                if (split.length != 2)
                    return paramError();

                gp.m_maxMissedCleavages = Integer.parseInt(split[1]);
            }
            else if (params[i].startsWith("-cys"))
            {
                String[] split = params[i].split("=");
                if (split.length != 2)
                    return paramError();
                double cys = Double.parseDouble(split[1]);
                double[] massTab = getMasses(true);
                massTab['C'] += cys;
                gp.m_massTab = massTab;
            }
            else if (!params[i].startsWith("-"))
            {
                gp.m_inputFileName = params[i];
            }
            else
            {
                //Unknown option
                return paramError();
            }
        }

        if (null == gp.m_inputFileName)
            return paramError();

        File file = new File(gp.m_inputFileName);
        if (!file.exists())
        {
            System.err.println("Could not find file: " + file.getAbsolutePath());
            return null;
        }

        gp.addListener(gp.getCommandLineListener());


        return gp;
    }

    private static PeptideGenerator paramError()
    {
        System.err.print("syntax: java PeptideGenerator inputFile [-compute={m,mm,pi,hp}] [-digest={tryptic|all}] [-massRange=min-max] [-countOnly] [-cys=m]");
        return null;
    }

    public void run()
    {
        FastaLoader loader = null;
        try
        {
            File fastaFile = new File(m_inputFileName);
            loader = new FastaLoader(fastaFile);
        }
        catch (Exception x)
        {
            x.printStackTrace(System.err);
            return;
        }

        Iterator it = loader.iterator();

        double[] massTab = getMasses(false);
        //it can be used for pruning purpose
        if (m_minResidues == -1)
            m_minResidues = (int) (m_minMass / massTab['W']);
        if (m_maxResidues == Integer.MAX_VALUE)
            m_maxResidues = (int) (m_maxMass / massTab['G']);

        while (it.hasNext())
        {
            m_protNum++;
            Protein p = (Protein) it.next();
            byte[] bytes = p.getBytes();
            m_aaCount += bytes.length;
            if (m_digest == DIGEST_TRYPTIC)
                doProteinDigest(p);
            else
                doProteinAll(p);

            if (m_protNum >= m_maxProteins)
                break;
        }

        fireHandleDone();

    }

    public Peptide[] digestProtein(Protein protein)
    {
        m_async = false;
        m_peptides = new ArrayList();

        doProteinDigest(protein);

        return (Peptide[]) m_peptides.toArray(new Peptide[m_peptides.size()]);
    }

    public double[] getPeptideMasses(Protein protein){
        
    	m_async = false;
    	m_peptides = new ArrayList();
        doProteinDigest(protein);
        //System.out.println(m_peptides.size());
        Peptide[] pep = (Peptide[]) m_peptides.toArray(new Peptide[m_peptides.size()]);
        
        /*Added by Zengyou. We consider the variable modification of Met-Oxidation: +16Da shift*/
        /*
        int[] numMet = new int[pep.length];
        int totalNumMet = 0;
        for(int i=0;i<pep.length;i++){
        	String str = pep[i].toString();
        	numMet[i] = 0;
        	for(int j=0;j<str.length();j++){
        		if(str.charAt(j)=='M'){
        			numMet[i]++;
        			totalNumMet++;
        		}
        	}
        }
        */
        //double[] masses = new double[pep.length+totalNumMet]; //changed by zengyou
        double[] masses = new double[pep.length];
        int counter = 0;
        for(int i=0;i<pep.length;i++){
        	/* We believe that "m_massTab" has been initialized to desired value*/
        	masses[i] = pep[i].getMass(m_massTab);
        	
        	/// added by zengyou
        	/*
        	int j=1;
        	while(j<=numMet[i]){
        		masses[pep.length+counter] = masses[i]+15.9945*j;
        		j++;
        		counter++;
        	}
        	*/
        	/// added by zengyou
        }
        return masses;
    }

    public double getProteinMass(Protein protein){
    	Peptide pep = new Peptide(protein,0,protein.getBytes().length-1);
    	return pep.getMass(m_massTab);
    }
    public void addListener(PeptideListener listener)
    {
        m_listeners.add(listener);
    }

    public PeptideListener getCommandLineListener()
    {
        return new CommandLinePeptideListener();
    }

    private void doProteinDigest(Protein protein)
    {
        byte[] bytes = protein.getBytes();
        ArrayList rkSites = new ArrayList();

        for (int i = 0; i < bytes.length; i++)
        {
            byte b = bytes[i];
            if ((i + 1 == bytes.length || bytes[i + 1] != 'P') && (b == 'R' || b == 'K'))
                    rkSites.add(new Integer(i));
        }
        byte b = bytes[bytes.length -1];
        if (b != 'K' && b != 'R')
            rkSites.add(new Integer(bytes.length - 1));

        Integer[] rk = (Integer[]) rkSites.toArray(new Integer[rkSites.size()]);
        int prevSite = -1;

        int nextSite = 0;
        for (int i = 0; i < rk.length; i++)
        {
            double m = 0;
            for (int j = 0; j <= m_maxMissedCleavages && i + j < rk.length; j++)
            {
                nextSite = rk[i + j].intValue();
                if (nextSite - prevSite > m_minResidues && nextSite - prevSite < m_maxResidues)
                {
                    Peptide pep = new Peptide(protein, prevSite + 1, nextSite - prevSite);
                    m = pep.getMass(m_massTab);
                    if (m >= m_minMass && m <= m_maxMass)
                        if (m_async)
                            fireHandlePeptide(pep);
                        else
                            m_peptides.add(pep);
                }
                if (m > m_maxMass)
                    break;
            }
            prevSite = rk[i].intValue();
        }
    }

    //Digestion is not performed  
    private void doProteinAll(Protein protein)
    {
        byte[] bytes = protein.getBytes();
        double[] massTab = getMasses(false);
        for (int i = 0; i < bytes.length - m_minResidues; i++)
        {
            double m = massTab['h'] + massTab['o'] + massTab['h'];

            for (int j = 0; j < m_minResidues; j++)
                m +=  massTab[bytes[i + j]];

            for (int j = m_minResidues; j <= m_maxResidues && i + j < bytes.length; j++)
            {
                m += massTab[bytes[i + j]];
                //m = fireHandlePeptide(bytes, i, j);
                if (m > m_maxMass)
                    break;
                else if (m > m_minMass)
                {
                    if (m_countOnly)
                        m_pepCount++;
                    else
                        fireHandlePeptide(new Peptide(protein, i, j + 1));
                }
                if(m_pepCount > 0 && m_pepCount % 1000000 == 0)
                    System.err.println(m_pepCount);
            }
        }
    }


    private void fireHandlePeptide(Peptide peptide)
    {
        Iterator iter = m_listeners.iterator();
        while (iter.hasNext())
        {
            PeptideListener listener = (PeptideListener) iter.next();
            listener.handlePeptide(peptide);
        }
    }

    private void fireHandleDone()
    {
        Iterator iter = m_listeners.iterator();
        while (iter.hasNext())
        {
            PeptideListener listener = (PeptideListener) iter.next();
            listener.handleDone();
        }
    }



    public static double computeMass(byte[] bytes, int start, int length, double[] massTab)
    {
        double pepMass = massTab['h'] + massTab['o'] + massTab['h'];
        for (int a = start; a < start + length; a++ )
            pepMass += massTab[bytes[a]];

        return pepMass;
    }

    
    /**
     * Taken from AminoAcidMasses.h on sourceforge.net. Returns 128 character array.
     *
     * Lower case indexes h, o, c, n, p, s are masses for corresponding elements
     * Upper case indexes are Amino Acid Masses.
     *
     * @param monoisotopic If true, return monoisotopic masses, otherwise average masses
     * @return
     */
    public static double[] getMasses(
            boolean monoisotopic)
    {
       double[] aaMasses = new double[128];
       if (!monoisotopic)
       {
          aaMasses['h']=  1.00794;  /* hydrogen */
          aaMasses['o']= 15.9994;   /* oxygen */
          aaMasses['c']= 12.0107;   /* carbon */
          aaMasses['n']= 14.00674;  /* nitrogen */
          aaMasses['p']= 30.973761; /* phosporus */
          aaMasses['s']= 32.066;    /* sulphur */

          aaMasses['G']= 57.05192;
          aaMasses['A']= 71.07880;
          aaMasses['S']= 87.07820;
          aaMasses['P']= 97.11668;
          aaMasses['V']= 99.13256;
          aaMasses['T']=101.10508;
          aaMasses['C']=103.13880;//+57.0513; /* 103.1448, 103.14080 */ //Carboxamidomethyl(on Cysteine);Added by Zengyou(28-02-2008)
          aaMasses['L']=113.15944;
          aaMasses['I']=113.15944;
          aaMasses['X']=113.15944;
          aaMasses['N']=114.10384;
          aaMasses['O']=114.14720;
          aaMasses['B']=114.59622;
          aaMasses['D']=115.08860;
          aaMasses['Q']=128.13072;
          aaMasses['K']=128.17408;
          aaMasses['Z']=128.62310;
          aaMasses['E']=129.11548;
          aaMasses['M']=131.19256; /* 131.19456 131.1986 */
          aaMasses['H']=137.14108;
          aaMasses['F']=147.17656;
          aaMasses['R']=156.18748;
          aaMasses['Y']=163.17596;
          aaMasses['W']=186.21320;
       }
       else /* monoisotopic masses */
       {
          aaMasses['h']=  1.0078250;
          aaMasses['o']= 15.9949146;
          aaMasses['c']= 12.0000000;
          aaMasses['n']= 14.0030740;
          aaMasses['p']= 30.9737633;
          aaMasses['s']= 31.9720718;

          aaMasses['G']= 57.0214636;
          aaMasses['A']= 71.0371136;
          aaMasses['S']= 87.0320282;
          aaMasses['P']= 97.0527636;
          aaMasses['V']= 99.0684136;
          aaMasses['T']=101.0476782;
          aaMasses['C']=103.0091854;//+57.021464; //Carboxamidomethyl (on Cysteine); Added by Zengyou (28-02-2008)
          aaMasses['L']=113.0840636;
          aaMasses['I']=113.0840636;
          aaMasses['X']=113.0840636;
          aaMasses['N']=114.0429272;
          aaMasses['O']=114.0793126;
          aaMasses['B']=114.5349350;
          aaMasses['D']=115.0269428;
          aaMasses['Q']=128.0585772;
          aaMasses['K']=128.0949626;
          aaMasses['Z']=128.5505850;
          aaMasses['E']=129.0425928;
          aaMasses['M']=131.0404854;
          aaMasses['H']=137.0589116;
          aaMasses['F']=147.0684136;
          aaMasses['R']=156.1011106;
          aaMasses['Y']=163.0633282;
          aaMasses['W']=186.0793126;
       }

        aaMasses[H_ION_INDEX] = aaMasses['h'] - ELECTRON_MASS;  //why???
        return aaMasses;
    } /*ASSIGN_MASS*/


    /**
     * PI Calculation
     */

    private static final double PH_MIN = 0;       /* minimum pH value */
    private static final double PH_MAX =14;      /* maximum pH value */
    private static final double MAXLOOP = 2000;    /* maximum number of iterations */
    private static final double  EPSI   = 0.0001;  /* desired precision */


    /* the 7 amino acid which matter */
        static int R = 'R' - 'A',
                   H = 'H' - 'A',
                   K = 'K' - 'A',
                   D = 'D' - 'A',
                   E = 'E' - 'A',
                   C = 'C' - 'A',
                   Y = 'Y' - 'A';

    /*
     *  table of pk values :
     *  Note: the current algorithm does not use the last two columns. Each
     *  row corresponds to an amino acid starting with Ala. J, O and U are
     *  inexistant, but here only in order to have the complete alphabet.
     *
     *          Ct    Nt   Sm     Sc     Sn
     */

        static double[][]  pk = new double[][]
        {
/* A */    {3.55, 7.59, 0.0  , 0.0  , 0.0   },
/* B */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* C */    {3.55, 7.50, 9.00 , 9.00 , 9.00  },
/* D */    {4.55, 7.50, 4.05 , 4.05 , 4.05  },
/* E */    {4.75, 7.70, 4.45 , 4.45 , 4.45  },
/* F */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* G */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* H */    {3.55, 7.50, 5.98 , 5.98 , 5.98  },
/* I */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* J */    {0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* K */    {3.55, 7.50, 10.00, 10.00, 10.00 },
/* L */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* M */    {3.55, 7.00, 0.0  , 0.0  , 0.0   },
/* N */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* O */    {0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* P */    {3.55, 8.36, 0.0  , 0.0  , 0.0   },
/* Q */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* R */    {3.55, 7.50, 12.0 , 12.0 , 12.0  },
/* S */    {3.55, 6.93, 0.0  , 0.0  , 0.0   },
/* T */    {3.55, 6.82, 0.0  , 0.0  , 0.0   },
/* U */    {0.00, 0.00, 0.0  , 0.0  , 0.0   },
/* V */    {3.55, 7.44, 0.0  , 0.0  , 0.0   },
/* W */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* X */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
/* Y */    {3.55, 7.50, 10.00, 10.00, 10.00 },
/* Z */    {3.55, 7.50, 0.0  , 0.0  , 0.0   },
        };

        static double exp10(double value)
        {
           return Math.pow(10.0,value);
        }

        static double computePI(byte[] seq, int start, int seq_length)
        {
           int[]             comp = new int[26];    /* Amino acid composition of the uni.tubingen.protein */
           int    nterm_res,   /* N-terminal residue */
                   cterm_res;   /* C-terminal residue */
           int i;
            int charge_increment = 0;
           double charge,
                           ph_mid = 0,
                           ph_min,
                           ph_max;
           double          cter,
                           nter;
           double          carg,
                           clys,
                           chis,
                           casp,
                           cglu,
                           ctyr,
                           ccys;


           for (i = 0; i < seq_length; i++)        /* compute the amino acid composition */
           {
              comp[seq[i + start] - 'A']++;
           }

           nterm_res = seq[start] - 'A';               /* Look up N-terminal residue */
           cterm_res = seq[start + seq_length-1] - 'A';    /* Look up C-terminal residue */

           ph_min = PH_MIN;
           ph_max = PH_MAX;

           for (i = 0, charge = 1.0; i<MAXLOOP && (ph_max - ph_min)>EPSI; i++)
           {
              ph_mid = ph_min + (ph_max - ph_min) / 2.0;

              cter = exp10(-pk[cterm_res][0]) / (exp10(-pk[cterm_res][0]) + exp10(-ph_mid));
              nter = exp10(-ph_mid) / (exp10(-pk[nterm_res][1]) + exp10(-ph_mid));

              carg = comp[R] * exp10(-ph_mid) / (exp10(-pk[R][2]) + exp10(-ph_mid));
              chis = comp[H] * exp10(-ph_mid) / (exp10(-pk[H][2]) + exp10(-ph_mid));
              clys = comp[K] * exp10(-ph_mid) / (exp10(-pk[K][2]) + exp10(-ph_mid));

              casp = comp[D] * exp10(-pk[D][2]) / (exp10(-pk[D][2]) + exp10(-ph_mid));
              cglu = comp[E] * exp10(-pk[E][2]) / (exp10(-pk[E][2]) + exp10(-ph_mid));

              ccys = comp[C] * exp10(-pk[C][2]) / (exp10(-pk[C][2]) + exp10(-ph_mid));
              ctyr = comp[Y] * exp10(-pk[Y][2]) / (exp10(-pk[Y][2]) + exp10(-ph_mid));

              charge = carg + clys + chis + nter + charge_increment
                 - (casp + cglu + ctyr + ccys + cter);

              if (charge > 0.0)
              {
                 ph_min = ph_mid;
              }
              else
              {
                 ph_max = ph_mid;
              }
           }

           return ph_mid;
        }


    public String getInputFileName()
    {
        return m_inputFileName;
    }

    public void setInputFileName(String inputFileName)
    {
        this.m_inputFileName = inputFileName;
    }

    public String getOutputFileName()
    {
        return m_outputFileName;
    }

    public void setOutputFileName(String outputFileName)
    {
        this.m_outputFileName = outputFileName;
    }

    public boolean isComputePI()
    {
        return m_computePI;
    }

    public void setComputePI(boolean computePI)
    {
        this.m_computePI = computePI;
    }

    public boolean isComputeHp()
    {
        return m_computeHp;
    }

    public void setComputeHp(boolean computeHp)
    {
        this.m_computeHp = computeHp;
    }

    public int getHpWindowSize()
    {
        return m_hpWindowSize;
    }

    public void setHpWindowSize(int hpWindowSize)
    {
        this.m_hpWindowSize = hpWindowSize;
    }

    public boolean isComputeAverageMass()
    {
        return m_computeAverageMass;
    }

    public void setComputeAverageMass(boolean computeAverageMass)
    {
        this.m_computeAverageMass = computeAverageMass;
    }

    public boolean isComputeMonoisotopicMass()
    {
        return m_computeMonoisotopicMass;
    }

    public void setComputeMonoisotopicMass(boolean computeMonoisotopicMass)
    {
        this.m_computeMonoisotopicMass = computeMonoisotopicMass;
    }

    public int getDigest()
    {
        return m_digest;
    }

    public void setDigest(int digest)
    {
        this.m_digest = digest;
    }

    public boolean isCountOnly()
    {
        return m_countOnly;
    }

    public void setCountOnly(boolean countOnly)
    {
        this.m_countOnly = countOnly;
    }

    public int getMaxMissedCleavages()
    {
        return m_maxMissedCleavages;
    }

    public void setMaxMissedCleavages(int maxMissedCleavages)
    {
        m_maxMissedCleavages = maxMissedCleavages;
    }

    public double getMinMass()
    {
        return m_minMass;
    }

    public void setMinMass(double minMass)
    {
        m_minMass = minMass;
    }

    public double getMaxMass()
    {
        return m_maxMass;
    }

    public void setMaxMass(double maxMass)
    {
        m_maxMass = maxMass;
    }

    public int getMinResidues()
    {
        return m_minResidues;
    }

    public void setMinResidues(int minResidues)
    {
        m_minResidues = minResidues;
    }

    public int getMaxResidues()
    {
        return m_maxResidues;
    }

    public void setMaxResidues(int maxResidues)
    {
        m_maxResidues = maxResidues;
    }

    public double[] getMassTable()
    {
        return m_massTab;
    }

    public void setMassTable(double[] massTab)
    {
        m_massTab = massTab;
    }

    public class CommandLinePeptideListener implements PeptideListener
    {
        PrintStream out;

        public CommandLinePeptideListener()
        {
            out = System.out;
        }

        public void handlePeptide(Peptide peptide)
        {
            double m = 0;

            m = peptide.getMass();
            if (m > m_minMass && m <= m_maxMass)
            {
                m_pepCount++;

                if (!m_countOnly)
                {
                    out.print(m_protNum);
                    out.print('\t');

                    out.print(peptide.getChars());

                    if (m_computeAverageMass)
                    {
                        out.print('\t');
                        out.print(m);
                    }
                    if (m_computeMonoisotopicMass)
                    {
                        out.print('\t');
                        out.print(peptide.getMonoisotopicMass());
                    }
                    if (m_computePI)
                    {
                        out.print('\t');
                        out.print(peptide.getPi());
                    }
                    if(m_computeHp)
                    {
                        out.print('\t');
                        out.print(peptide.getHydrophobicity());
                    }
                    out.println();
                }
            }
        }

        public void handleDone()
        {

        }
    }

    public static interface PeptideListener
    {
        public void handlePeptide(Peptide peptide);
        public void handleDone();
    }
}
