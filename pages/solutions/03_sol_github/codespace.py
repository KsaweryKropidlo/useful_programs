def oblicz_srednia(wartosci):
	if not wartosci:
		return None
	return sum(wartosci) / len(wartosci)


def main():
	temperatury = [18.5, 21.0, 19.5, 22.0, 20.0]
	srednia = oblicz_srednia(temperatury)

	print("Pomiary temperatury:", temperatury)
	if srednia is None:
		print("Brak pomiarow do obliczenia sredniej.")
		return

	print(f"Srednia temperatura: {srednia:.1f} C")
	print(f"Najwyzsza temperatura: {max(temperatury):.1f} C")


if __name__ == "__main__":
	main()
