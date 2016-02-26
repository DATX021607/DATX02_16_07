import java.io.BufferedReader;
import java.io.FileReader;

public class Main {
	static BufferedReader br;
	public static void main(String args[]) {
		for(int i = 0; i < args.length; i++) {
			System.out.println(args[i].substring(0, 4));
			try {
				br = new BufferedReader(new FileReader(args[i]));
			    String line = br.readLine();
			    new FrequenceAnalysis(line);
			    br.close();
			} catch (Exception e) {
				
			}
		}
	}
}
