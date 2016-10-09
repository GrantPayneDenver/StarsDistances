/*
Reads in information from the HYG Database

The HYG (Hipparcos, Yale, Gliese) Database (v2.0) is a compilation of interesting stellar data from a
variety of catalogs. It is useful for background information on all sorts of data: star names, positions,
brightnesses, distances, and spectrum information. The HYG Database web site is at
http://www.astronexus.com/node/34.
The current version of the HYG Database is hosted at Github – https://github.com/astronexus/HYGDatabase.

The hygxyz.csv file is a text file in comma separated value (CSV) format. The first line of the file has the
names of the fields – 23 in all. These fields are: StarID, HIP, HD, HR, Gliese, BayerFlamsteed,
ProperName, RA, Dec, Distance, PMRA, PMDec, RV, Mag, AbsMag, Spectrum, ColorIndex, X, Y, Z, VX, VY,
and VZ. The remaining 119617 lines contain the stellar data – one star per line

Code reads in the X Y Z coordinates and computes information on them for each star.

This is a Java program to read the stars from the HYG database, ignoring those whose distances are not
accurately known, and compute the minimum, maximum, and mean shortest distances between stars.
That is, for each star, find the distance to its nearest neighbor and tally the statistics for that distance to
its nearest neighbor for all stars.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


public class Stars {

    public static double[][] starMatrix = new double[119617][3]; // stars by x y z coors

    // class holds the information that each starsifter gathers
    // once all are done, it prints out their information in the final format.
    public static class dataContainer
    {
        public double[] minAverages = new double [4];
        public double[] maxAdjacents = new double[4];
        public double[] minAdjacents = new double[4];
        int i =0;
        // method to input data

        synchronized void inputData(double minAverage, double maxAdjacent, double min, int id){
            minAverages[i] = minAverage;
            maxAdjacents[i] = maxAdjacent;
            minAdjacents[i] = min;
            ++i;
            if(i ==4)
            {
                double average = 0;
                double actualmax=0; double actualmin=1000000;
                for(double t: minAverages){average+=t;}
                for(double t: maxAdjacents){if(t > actualmax) actualmax = t;}
                for(double t: minAdjacents){if(t < actualmin) actualmin = t;}

                System.out.println("Average minmum distance: " + average/4);
                System.out.println("Minimum = : "+ actualmin);
                System.out.println("Maximum = : "+ actualmax);
            }

        }

    }

    public static void main(String[]args)
    {

        String line;
        BufferedReader br;
        int starCount=0;                    //sc = n - 1 elements

        try {
            br = new BufferedReader(new FileReader("C:\\Users\\Folio\\Desktop\\hygxyz.csv"));

            String throwAwayLine = br.readLine();
            String[] splitArray;

        while ((line = br.readLine()) !=null )     // reading in all stars
        {

            splitArray = line.split(",");

            if((!splitArray[9].equals("10000000")) /*&& Double.parseDouble(splitArray[9])<=10 */) // if Distance is known, add star to star matrix with x y z data
            {
                for (int j = 0; j < 3; ++j)
                {
                    starMatrix[starCount][j]=Double.parseDouble(splitArray[j + 17]);
                }
                ++starCount;
            }

        } // while

        }catch(IOException e){System.out.print(e);}

        /////////////////////////////////////////////////////////////////////
        // FINISHED GRABBING DATA
        /////////////////////////////////////////////////////////////////////

        /*

        Divide the array into fourths
        Make 4 StarSifter threads, have each run through doing the same calculation
        Once one is finished have it input it's data into an array that's synchronized and then go offline

        p = 0, r = n-1
        b = floor(r/4)


        A1 gets Array[p..b-1]
        A2 gets Array[b..2b-1]
        A3 gets Array[2b..3b-1]
        A4 gets Array[3b..r]

        A1 = [b][3] A2 = [b][3] A3 = [b][3] A4 = [n - 3b][3]
         */
///*
        System.out.println(starCount + " is the number of stars read in."); // 115372

        int p = 0;
        int r = starCount;
        int b = Math.floorDiv(r,4);

        dataContainer dc = new dataContainer();

        Thread sift1 = new StarSifter(1,  p,    b-1,   r+1, dc); // params: ID num, array seg start, array seg end, starCount
        Thread sift2 = new StarSifter(2,  b,    2*b-1, r+1, dc);
        Thread sift3 = new StarSifter(3,  2*b,  3*b-1, r+1, dc);
        Thread sift4 = new StarSifter(4,  3*b,  r,     r+1, dc);

        System.out.println("b1 "+p+" b2 "+ (b-1));
        System.out.println("b1 "+b+" b2 "+((2*b)-1));
        System.out.println("b1 "+2*b+" b2 "+(3*b-1));
        System.out.println("b1 "+3*b+" b2 "+r);

        sift1.start();
        sift2.start();
        sift3.start();
        sift4.start();

    } // main

    public static class StarSifter extends Thread
    {
        private int bound1, bound2;
        private int starCount;
        private int ID;
        private int starsRun=0;
        boolean notdone = true;

        dataContainer dc;

        double m_cartesianDistance = 0;
        double m_shortestDist = 10000000;
        double m_thisStarMin;
        double m_largestMostAdjacent=0;
        double m_minAvg=0;


        public StarSifter(int ID, int startOfArraySegment,int endOfArraySegment, int starCount, dataContainer dc)
        {
            this.dc = dc;
            this.ID = ID;
            this.starCount = starCount;
            this.bound1 = startOfArraySegment;
            this.bound2 = endOfArraySegment;
        }

        public void run()
        {
            for(int star = bound1; star < bound2; ++star)                         // for each star
            {
                ++starsRun;
                m_thisStarMin = 100000000;                                        // reset these vals per star

                double x1 = starMatrix[star][0];
                double y1 = starMatrix[star][1];
                double z1 = starMatrix[star][2];

                for(int otherStar=0; otherStar<starCount; ++otherStar)            // for each other star
                {
                    if(star == otherStar) continue;                               // not comparing star to itself

                    double x2 = starMatrix[otherStar][0];
                    double y2 = starMatrix[otherStar][1];
                    double z2 = starMatrix[otherStar][2];

                    m_cartesianDistance = Math.sqrt(Math.pow((x1 - x2),2)+Math.pow((y1 - y2),2)+Math.pow((z1 - z2),2));

                    if(m_cartesianDistance < m_thisStarMin) m_thisStarMin = m_cartesianDistance;  // finding this star's closest neighbor
                    if(m_cartesianDistance < m_shortestDist) m_shortestDist = m_cartesianDistance;// absolute shortest distance


                } // for otherStar

                m_minAvg+=m_thisStarMin;                                                          // adding smallest dist to smallest dist average
                if(m_largestMostAdjacent < m_thisStarMin) m_largestMostAdjacent = m_thisStarMin;  // starts at zero, gets largest thisStarMin

            } // for star

            m_minAvg /= starsRun;

            notdone = false;

            dc.inputData(m_minAvg, m_largestMostAdjacent, m_shortestDist, ID);

        } // run()

    } // starSifter

} // Stars
