import logging
from rdkit import Chem
from rdkit.Chem import rdChemReactions

logger = logging.getLogger(__name__)


def get_ec_class_scenarios(ec_class="3.1.1.1"):
    logger.info(f"Запрос сценариев для EC-класса: {ec_class}")
    if ec_class.startswith("3.1"):
        return {
            "mechanism": "Ser-His-Asp Triad",
            "pharmacophores": ["Ser_OG", "His_NE2", "Asp_OD1", "Oxyanion_Hole"],
            "target_bond": "ester",
        }
    logger.warning(
        f"EC-класс {ec_class} не опознан, используется стандартный сценарий."
    )
    return {"mechanism": "Unknown", "pharmacophores": [], "target_bond": "any"}


def get_retrosynthetic_network(smiles, ec_class="3.1.1.1"):
    logger.info(f"Запуск ретросинтеза для SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error("Не удалось распарсить SMILES!")
        return []

    scenarios = get_ec_class_scenarios(ec_class)
    intermediates = []

    if scenarios["target_bond"] == "ester":
        # Реакция гидролиза сложного эфира
        rxn_smarts = (
            "[CX3:1](=[OX1:2])[OX2:3][#6:4]>>[CX3:1](=[OX1:2])[OH].[OX2H:3][#6:4]"
        )
        logger.debug(f"Используем SMARTS для гидролиза: {rxn_smarts}")

        try:
            rxn = rdChemReactions.ReactionFromSmarts(rxn_smarts)
            products_sets = rxn.RunReactants((mol,))
            logger.info(f"Найдено возможных путей расщепления: {len(products_sets)}")

            for i, prods in enumerate(products_sets):
                if len(prods) == 2:
                    prod1_smiles = Chem.MolToSmiles(prods[0])
                    prod2_smiles = Chem.MolToSmiles(prods[1])
                    final_prod = f"{prod1_smiles}.{prod2_smiles}"

                    intermediates.append(
                        {
                            "path_id": i + 1,
                            "site": f"Связь {i + 1} (sn-{i + 1})",
                            "intermediate_type": "Тетраэдрический интермедиат",
                            "final_products": final_prod,
                            "activation_energy_shift": i * 5.0,
                        }
                    )
                    logger.debug(f"Путь {i + 1} успешно обработан: {final_prod}")
                else:
                    logger.warning(
                        f"Путь {i + 1} вернул неожиданное количество продуктов: {len(prods)}"
                    )

        except Exception as e:
            logger.error(f"Ошибка RDKit при ретросинтезе: {e}", exc_info=True)

    logger.info(f"Ретросинтез завершен. Возвращено путей: {len(intermediates)}")
    return intermediates
