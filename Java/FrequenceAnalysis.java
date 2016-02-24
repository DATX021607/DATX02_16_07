import java.util.HashMap;
import java.util.Map;

public class FrequenceAnalysis {
	private  Map<Character, Integer> map = new HashMap<Character,Integer>();
	public FrequenceAnalysis(String string) {
		String input = string;
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
	}
}
