import java.util.HashMap;
import java.util.Map;
import java.util.Scanner;

public class FrequenceAnalysis {

	private Map<Character, Integer> map = new HashMap<Character,Integer>();
	public FrequenceAnalysis() {
		Scanner sc = new Scanner(System.in);
		System.out.println("Input the sequence:");
		String input = sc.nextLine();
		char[] characters = input.toCharArray();
		for(char c : characters) {
			if(map.containsKey(c)) {
				map.replace(c, map.get(c) + 1);
			} else {
				map.put(c, 1);
			}
		}
		System.out.println(map.keySet());
		System.out.println(map.values());
		sc.close();
	}
}
